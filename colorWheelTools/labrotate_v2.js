/* =====================================================================
 * labrotate_wheel.js — object colour rotation onto a fixed CIELAB wheel
 *
 * Differences from labrotate.js (the RotateExplicitPar.m port):
 *
 *  1. The target is an explicit wheel  W(theta) = (L0, a0 + r cos theta,
 *     b0 + r sin theta)  — Zhang & Luck (2008) style — instead of "rotate
 *     every pixel's a*b* about the neutral axis". theta is now an ABSOLUTE
 *     wheel position, the same number you use for the response wheel.
 *
 *  2. Each object is anchored to the wheel by its mean colour:
 *         lab'(p) = W(theta) + R(theta - h0) * (lab(p) - mean)
 *     The object's mean a*b* lands exactly on W(theta); each pixel's
 *     deviation from the mean (shading, texture, highlights) is carried
 *     along and rotated with it. Optionally the mean L* is moved to L0.
 *
 *  3. Correct sRGB transfer on BOTH ends. The original decodes with a
 *     Rec.709-style curve and does not re-encode on output
 *     (APPLY_OUTPUT_GAMMA = false), so linear values went straight to the
 *     framebuffer and the on-screen colours were not the Lab values you
 *     computed. Set LEGACY_TRANSFER = true to reproduce the old behaviour.
 *
 *  4. Out-of-gamut pixels are counted per frame (this.clipFrac) so you can
 *     check how often the gamut mapper is altering colours.
 * ===================================================================== */

'use strict';

/* ---- wheel presets --------------------------------------------------- */
const WHEELS = {
  // Zhang & Luck (2008). Calibrated CRT; about 200 of 360 degrees fall
  // OUTSIDE the sRGB gamut, so it cannot be shown faithfully on a
  // standard / uncalibrated display.
  zhangLuck: { L: 70, a0: 20,   b0: 38,   r: 60 },
  // Schurgin, Wixted & Brady (2020) and later Brady-lab studies.
  // Entirely inside sRGB. Recommended for online / uncalibrated work.
  schurgin:  { L: 54, a0: 21.5, b0: 11.5, r: 49 },
};

const LEGACY_TRANSFER = false;   // true = old Rec.709 in / no gamma out

/* ---- colorimetry (unchanged from the MATLAB port) -------------------- */
const WP0 = 0.950456, WP1 = 1, WP2 = 1.088754;   // D65, 2-deg

const M00 =  3.240479, M01 = -1.537150, M02 = -0.498535;
const M10 = -0.969256, M11 =  1.875992, M12 =  0.041556;
const M20 =  0.055648, M21 = -0.204043, M22 =  1.057311;

const Mi00 = 0.412453, Mi01 = 0.357580, Mi02 = 0.180423;
const Mi10 = 0.212671, Mi11 = 0.715160, Mi12 = 0.072169;
const Mi20 = 0.019334, Mi21 = 0.119193, Mi22 = 0.950227;

const K1 = 841 / 108, K2 = 4 / 29, K3 = 108 / 841;

/* ---- transfer functions ---------------------------------------------- */
function legacyDecode(Rp) {                       // the .m file's curve
  let R = Math.pow((Rp + 0.099) / 1.099, 1 / 0.45);
  if (R < 0.018) R = Rp / 4.5138;
  return R;
}
function srgbDecode(v) {                          // IEC 61966-2-1
  return v <= 0.04045 ? v / 12.92 : Math.pow((v + 0.055) / 1.055, 2.4);
}
const decode = LEGACY_TRANSFER ? legacyDecode : srgbDecode;

// Linear -> 8-bit via a LUT so the hot loop stays transcendental-free.
const ENC_N = 4096;
const ENC_LUT = new Uint8ClampedArray(ENC_N + 1);
for (let i = 0; i <= ENC_N; i++) {
  const x = i / ENC_N;
  const v = LEGACY_TRANSFER ? x
          : (x <= 0.0031308 ? 12.92 * x : 1.055 * Math.pow(x, 1 / 2.4) - 0.055);
  ENC_LUT[i] = Math.round(255 * v);
}
function encode8(x) { return ENC_LUT[(x * ENC_N + 0.5) | 0]; }

function labf(Y)   { return Y < 0.008856 ? Y * K1 + K2 : Math.pow(Y, 1 / 3); }
function labinvf(f){ const Y = f * f * f; return Y < 0.008856 ? (f - K2) * K3 : Y; }

const D2R = Math.PI / 180;

/* Lab -> linear RGB with the original "add white, rescale" gamut mapping.
 * Writes into tmp[0..2]; returns true if the colour was out of gamut.   */
function labToLinear(L, A, B, tmp) {
  const fY = (L + 16) / 116;
  const X = WP0 * labinvf(fY + A / 500);
  const Y = WP1 * labinvf(fY);
  const Z = WP2 * labinvf(fY - B / 200);

  let R  = M00 * X + M01 * Y + M02 * Z;
  let G  = M10 * X + M11 * Y + M12 * Z;
  let Bc = M20 * X + M21 * Y + M22 * Z;

  let mn = R < G ? R : G; if (Bc < mn) mn = Bc;
  let mx = R > G ? R : G; if (Bc > mx) mx = Bc;
  const oog = mn < 0 || mx > 1;

  const addWhite = mn < 0 ? -mn : 0;
  let scale = mx + addWhite; if (scale < 1) scale = 1;
  const inv = 1 / scale;
  R = (R + addWhite) * inv; G = (G + addWhite) * inv; Bc = (Bc + addWhite) * inv;

  tmp[0] = R < 0 ? 0 : R > 1 ? 1 : R;
  tmp[1] = G < 0 ? 0 : G > 1 ? 1 : G;
  tmp[2] = Bc < 0 ? 0 : Bc > 1 ? 1 : Bc;
  return oog;
}

/* Colour of the wheel itself at angle theta (for painting the response
 * wheel and for any colour-patch conditions).                          */
function wheelRgb255(theta, wheel = WHEELS.schurgin) {
  const t = new Float64Array(3);
  labToLinear(wheel.L,
              wheel.a0 + wheel.r * Math.cos(theta * D2R),
              wheel.b0 + wheel.r * Math.sin(theta * D2R), t);
  return [encode8(t[0]), encode8(t[1]), encode8(t[2])];
}

/* ------------------------------------------------------------------- */
class LabStim {
  /**
   * @param {Uint8ClampedArray} rgba  straight-alpha RGBA
   * @param {number} w
   * @param {number} h
   * @param {object} [opts]
   * @param {object} [opts.wheel=WHEELS.schurgin]
   * @param {boolean} [opts.matchLightness=false]  shift mean L* to wheel.L
   * @param {number}  [opts.devScale=1]  scale deviations from the mean
   *        (<1 reduces within-object colour spread, fewer gamut clips)
   */
  constructor(rgba, w, h, opts = {}) {
    this.wheel = opts.wheel || WHEELS.schurgin;
    this.matchLightness = !!opts.matchLightness;
    this.devScale = opts.devScale ?? 1;

    this.width = w; this.height = h;
    const n = w * h; this.n = n;
    this.L = new Float64Array(n);
    this.A = new Float64Array(n);
    this.B = new Float64Array(n);
    this.alpha = new Uint8ClampedArray(n);
    this.out = new Uint8ClampedArray(n * 4);

    const idx = [];
    let sL = 0, sA = 0, sB = 0, sW = 0;

    for (let p = 0; p < n; p++) {
      const o = p << 2, a8 = rgba[o + 3];
      this.alpha[p] = a8;
      if (a8 === 0) continue;
      idx.push(p);

      const a = a8 / 255;
      const R  = decode((rgba[o]     / 255) * a + 0.5 * (1 - a));
      const G  = decode((rgba[o + 1] / 255) * a + 0.5 * (1 - a));
      const Bl = decode((rgba[o + 2] / 255) * a + 0.5 * (1 - a));

      const X = Mi00 * R + Mi01 * G + Mi02 * Bl;
      const Y = Mi10 * R + Mi11 * G + Mi12 * Bl;
      const Z = Mi20 * R + Mi21 * G + Mi22 * Bl;

      const fX = labf(X / WP0), fY = labf(Y / WP1), fZ = labf(Z / WP2);
      const L = 116 * fY - 16, A = 500 * (fX - fY), B = 200 * (fY - fZ);
      this.L[p] = L; this.A[p] = A; this.B[p] = B;

      // alpha-weighted mean: semi-transparent edge pixels are mostly grey
      sL += a * L; sA += a * A; sB += a * B; sW += a;
    }
    this.idx = Int32Array.from(idx);
    for (let p = 0; p < n; p++) this.out[(p << 2) + 3] = this.alpha[p];

    sW = sW || 1;
    this.meanL = sL / sW; this.meanA = sA / sW; this.meanB = sB / sW;
    // hue of the object's own mean colour; deviations are rotated by
    // (theta - h0) so the texture keeps its original relation to the mean
    this.h0 = Math.atan2(this.meanB, this.meanA) / D2R;

    this.lastTheta = null;
    this.clipFrac = 0;
    this._t = new Float64Array(3);
  }

  /** Render the object at absolute wheel angle theta (degrees). */
  render(theta) {
    if (theta === this.lastTheta) return this.out;
    this.lastTheta = theta;

    const W = this.wheel;
    const ta = W.a0 + W.r * Math.cos(theta * D2R);
    const tb = W.b0 + W.r * Math.sin(theta * D2R);
    const dL = this.matchLightness ? W.L - this.meanL : 0;

    const k = this.devScale;
    const phi = (theta - this.h0) * D2R;
    const c = k * Math.cos(phi), s = k * Math.sin(phi);

    const { L, A, B, idx, out, meanA, meanB } = this;
    const t = this._t;
    let clipped = 0;

    for (let i = 0; i < idx.length; i++) {
      const p = idx[i];
      const da = A[p] - meanA, db = B[p] - meanB;
      const Av = ta + c * da - s * db;
      const Bv = tb + s * da + c * db;

      if (labToLinear(L[p] + dL, Av, Bv, t)) clipped++;

      const o = p << 2;
      out[o]     = encode8(t[0]);
      out[o + 1] = encode8(t[1]);
      out[o + 2] = encode8(t[2]);
    }
    this.clipFrac = idx.length ? clipped / idx.length : 0;
    return out;
  }
}

if (typeof module !== 'undefined')
  module.exports = { LabStim, WHEELS, wheelRgb255, labToLinear };
