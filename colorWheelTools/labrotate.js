/* =====================================================================
 * labrotate.js — browser port of RotateExplicitPar.m
 *
 * Reproduces the MATLAB pipeline exactly:
 *   base = rgb.*alpha + 0.5.*(1-alpha)      composite onto neutral grey
 *   lab  = rgb2lab_explicit(base)           ONCE per stimulus
 *   rot  = [L, c*a - s*b, s*a + c*b]        per requested angle
 *   out  = lab2rgb_explicit(rot, false)     APPLY_OUTPUT_GAMMA = false
 *
 * The inverse path contains no transcendental functions (gamma is
 * disabled), so the per-frame cost is plain arithmetic on typed arrays.
 * ===================================================================== */

'use strict';

const WP0 = 0.950456, WP1 = 1, WP2 = 1.088754;   // D65, 2-deg observer

// XYZ -> linear sRGB (M in the .m file)
const M00 =  3.240479, M01 = -1.537150, M02 = -0.498535;
const M10 = -0.969256, M11 =  1.875992, M12 =  0.041556;
const M20 =  0.055648, M21 = -0.204043, M22 =  1.057311;

// linear sRGB -> XYZ.  MATLAB computes this numerically as inv(M).
// Run  `format long; inv(M)`  and paste your own digits if you want to be
// certain; these agree with inv(M) to ~1e-7.
const Mi00 = 0.412453, Mi01 = 0.357580, Mi02 = 0.180423;
const Mi10 = 0.212671, Mi11 = 0.715160, Mi12 = 0.072169;
const Mi20 = 0.019334, Mi21 = 0.119193, Mi22 = 0.950227;

const K1 = 841 / 108, K2 = 4 / 29, K3 = 108 / 841;

// MATLAB: R = real(((Rp+0.099)/1.099).^(1/0.45)); i = R < 0.018; R(i) = Rp(i)/4.5138
// NOTE the threshold is tested on the power-law RESULT, not on Rp.
function invgammacorrection(Rp) {
  let R = Math.pow((Rp + 0.099) / 1.099, 1 / 0.45);
  if (R < 0.018) R = Rp / 4.5138;
  return R;
}

// MATLAB: fY = real(Y.^(1/3)); i = Y < 0.008856; fY(i) = Y(i)*(841/108) + 4/29
// Threshold tested on Y, not on fY. Math.pow (not cbrt) to mirror Y.^(1/3).
function labf(Y) {
  if (Y < 0.008856) return Y * K1 + K2;
  return Math.pow(Y, 1 / 3);
}

// MATLAB: Y = fY.^3; i = Y < 0.008856; Y(i) = (fY(i) - 4/29)*(108/841)
// Threshold tested on the CUBE, not on fY. Handles negatives via the
// linear branch, so no NaN even when rotation pushes fX/fZ below zero.
function labinvf(f) {
  const Y = f * f * f;
  if (Y < 0.008856) return (f - K2) * K3;
  return Y;
}

// MATLAB cosd/sind are exact at multiples of 90; Math.cos(PI/2) is 6.1e-17.
// Irrelevant at 8-bit output but free to match.
function cosd(d) {
  const m = ((d % 360) + 360) % 360;
  if (m === 0) return 1;
  if (m === 90 || m === 270) return 0;
  if (m === 180) return -1;
  return Math.cos(m * Math.PI / 180);
}
function sind(d) {
  const m = ((d % 360) + 360) % 360;
  if (m === 0 || m === 180) return 0;
  if (m === 90) return 1;
  if (m === 270) return -1;
  return Math.sin(m * Math.PI / 180);
}

/* ---------------------------------------------------------------------
 * A single stimulus, converted to Lab once and then rotated on demand.
 * ------------------------------------------------------------------- */
class LabStim {
  /**
   * @param {Uint8ClampedArray} rgba  straight-alpha RGBA, length 4*w*h
   * @param {number} w
   * @param {number} h
   */
  constructor(rgba, w, h) {
    this.width = w;
    this.height = h;
    const n = w * h;
    this.n = n;

    this.L = new Float64Array(n);
    this.A = new Float64Array(n);
    this.B = new Float64Array(n);
    this.alpha = new Uint8ClampedArray(n);
    this.out = new Uint8ClampedArray(n * 4);

    // Pixels with alpha == 0 are invisible: their colour never reaches the
    // screen, so skip them in the per-frame loop. Typically 40-70% of a
    // cropped object bitmap.
    const idx = [];

    for (let p = 0; p < n; p++) {
      const o = p << 2;
      const a8 = rgba[o + 3];
      this.alpha[p] = a8;
      if (a8 === 0) continue;
      idx.push(p);

      const a = a8 / 255;
      // base = rgb.*a + 0.5.*(1-a)
      const br = (rgba[o]     / 255) * a + 0.5 * (1 - a);
      const bg = (rgba[o + 1] / 255) * a + 0.5 * (1 - a);
      const bb = (rgba[o + 2] / 255) * a + 0.5 * (1 - a);

      const R = invgammacorrection(br);
      const G = invgammacorrection(bg);
      const Bl = invgammacorrection(bb);

      const X = Mi00 * R + Mi01 * G + Mi02 * Bl;
      const Y = Mi10 * R + Mi11 * G + Mi12 * Bl;
      const Z = Mi20 * R + Mi21 * G + Mi22 * Bl;

      const fX = labf(X / WP0);
      const fY = labf(Y / WP1);
      const fZ = labf(Z / WP2);

      this.L[p] = 116 * fY - 16;
      this.A[p] = 500 * (fX - fY);
      this.B[p] = 200 * (fY - fZ);
    }

    this.idx = Int32Array.from(idx);
    for (let p = 0; p < n; p++) this.out[(p << 2) + 3] = this.alpha[p];

    this.lastTheta = null;
  }

  /**
   * Area-weighted mean lightness, chroma and hue angle over the visible
   * pixels. Lets the response wheel be built from this object's own colour
   * statistics instead of an arbitrary fixed reference.
   * @returns {{L:number, C:number, h:number}} h in degrees
   */
  meanChromaPolar() {
    const idx = this.idx, L = this.L, A = this.A, B = this.B;
    let sumL = 0, sumC = 0, sx = 0, sy = 0;
    for (let k = 0; k < idx.length; k++) {
      const p = idx[k];
      const a = A[p], b = B[p];
      const c = Math.hypot(a, b);
      sumL += L[p];
      sumC += c;
      sx += a;              // chroma-weighted already, since a = c*cos(h)
      sy += b;
    }
    const n = idx.length || 1;
    return {
      L: sumL / n,
      C: sumC / n,
      h: (Math.atan2(sy, sx) * 180 / Math.PI + 360) % 360
    };
  }

  /**
   * Rotate the a*b* plane by theta degrees and write RGBA into this.out.
   * @returns {Uint8ClampedArray} the RGBA buffer (reused, do not retain)
   */
  render(theta) {
    if (theta === this.lastTheta) return this.out;
    this.lastTheta = theta;

    const c = cosd(theta), s = sind(theta);
    const L = this.L, A = this.A, B = this.B;
    const idx = this.idx, out = this.out;
    const nIdx = idx.length;

    for (let k = 0; k < nIdx; k++) {
      const p = idx[k];
      const a0 = A[p], b0 = B[p];

      const Av = c * a0 - s * b0;
      const Bv = s * a0 + c * b0;

      const fY = (L[p] + 16) / 116;
      const fX = fY + Av / 500;
      const fZ = fY - Bv / 200;

      const X = WP0 * labinvf(fX);
      const Y = WP1 * labinvf(fY);
      const Z = WP2 * labinvf(fZ);

      let R  = M00 * X + M01 * Y + M02 * Z;
      let G  = M10 * X + M11 * Y + M12 * Z;
      let Bc = M20 * X + M21 * Y + M22 * Z;

      // out-of-gamut: add white to lift the negative channel, then rescale
      let mn = R < G ? R : G; if (Bc < mn) mn = Bc;
      const addWhite = mn < 0 ? -mn : 0;

      let mx = R > G ? R : G; if (Bc > mx) mx = Bc;
      let scale = mx + addWhite; if (scale < 1) scale = 1;
      const inv = 1 / scale;

      R  = (R  + addWhite) * inv;
      G  = (G  + addWhite) * inv;
      Bc = (Bc + addWhite) * inv;

      if (R < 0) R = 0; else if (R > 1) R = 1;
      if (G < 0) G = 0; else if (G > 1) G = 1;
      if (Bc < 0) Bc = 0; else if (Bc > 1) Bc = 1;

      const o = p << 2;
      // Math.round matches MATLAB's im2uint8 (half away from zero);
      // bare assignment to Uint8ClampedArray would round half to even.
      out[o]     = Math.round(255 * R);
      out[o + 1] = Math.round(255 * G);
      out[o + 2] = Math.round(255 * Bc);
    }
    return out;
  }
}


/* ---------------------------------------------------------------------
 * Standalone Lab -> sRGB, identical math to LabStim.render's inner loop.
 * Kept separate (rather than factored out) so the hot loop stays inlined.
 * Used to paint the response wheel in the same colour space as the stimuli,
 * including the disabled output gamma.
 * ------------------------------------------------------------------- */
function labToRgb255(L, A, B) {
  const fY = (L + 16) / 116;
  const fX = fY + A / 500;
  const fZ = fY - B / 200;

  const X = WP0 * labinvf(fX);
  const Y = WP1 * labinvf(fY);
  const Z = WP2 * labinvf(fZ);

  let R  = M00 * X + M01 * Y + M02 * Z;
  let G  = M10 * X + M11 * Y + M12 * Z;
  let Bc = M20 * X + M21 * Y + M22 * Z;

  let mn = R < G ? R : G; if (Bc < mn) mn = Bc;
  const addWhite = mn < 0 ? -mn : 0;
  let mx = R > G ? R : G; if (Bc > mx) mx = Bc;
  let scale = mx + addWhite; if (scale < 1) scale = 1;
  const inv = 1 / scale;

  R  = (R  + addWhite) * inv;
  G  = (G  + addWhite) * inv;
  Bc = (Bc + addWhite) * inv;

  if (R < 0) R = 0; else if (R > 1) R = 1;
  if (G < 0) G = 0; else if (G > 1) G = 1;
  if (Bc < 0) Bc = 0; else if (Bc > 1) Bc = 1;

  return [Math.round(255 * R), Math.round(255 * G), Math.round(255 * Bc)];
}

if (typeof module !== 'undefined') module.exports = { LabStim, labToRgb255, invgammacorrection, labf, labinvf };
