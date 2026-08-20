function RotateExplicitPar()
% ROTATEEXPLICITPAR  Parallel version of RotateExplicit.
%
%   1. parfor over images (coarse-grained, one worker per stimulus).
%   2. unique-colour LUT: the Lab->RGB round trip runs once per DISTINCT
%      colour per angle instead of once per pixel per angle.
%   3. outputs are named by each stimulus's ORIGINAL file ID, not by its
%      position in the dir() listing.
%
% Numerics are bit-identical to RotateExplicit (same constants, same
% quirky threshold tests, same disabled output gamma).

    APPLY_OUTPUT_GAMMA = false;   % false = match colorspace.m exactly
    USE_COLOR_LUT      = true;    % set false to verify against serial output
    ID_MODE            = "name";  % "name"   -> use full original basename
                                  % "digits" -> pull first digit run, zero-pad
    THETAS             = 1:359;

    stimSet = '/Users/ali/Desktop/MyExperiment/TestObjectsTransparent';
    outRoot = '/Users/ali/Desktop/MyExperiment/studyObjectsRotatedVars';

    stimFiles = dir(fullfile(stimSet, '*.png'));
    stimFiles = stimFiles(~[stimFiles.isdir]);
    stimNames = string({stimFiles.name})';
    stimPaths = fullfile(string(stimSet), stimNames);

    % --- derive a stable ID per stimulus from its original filename -------
    [~, baseNames] = fileparts(stimNames);

    switch ID_MODE
        case "name"
            stimIDs = baseNames;
        case "digits"
            nums = str2double(regexp(baseNames, '\d+', 'match', 'once'));
            if any(isnan(nums))
                bad = baseNames(isnan(nums));
                error('RotateExplicitPar:noDigits', ...
                      'No digits found in: %s', strjoin(bad, ', '));
            end
            stimIDs = compose("%03d", nums);
        otherwise
            error('RotateExplicitPar:badIDMode', 'Unknown ID_MODE.');
    end

    % IDs become directory names, so collisions would silently merge output.
    if numel(unique(stimIDs)) ~= numel(stimIDs)
        [u, ~, g] = unique(stimIDs);
        dupIDs = u(accumarray(g, 1) > 1);
        error('RotateExplicitPar:dupID', ...
              'Duplicate stimulus IDs: %s', strjoin(dupIDs, ', '));
    end

    % --- write a manifest so ID -> source file stays traceable ------------
    if ~isfolder(outRoot); mkdir(outRoot); end
    writetable(table(stimIDs, stimNames, 'VariableNames', {'id', 'sourceFile'}), ...
               fullfile(outRoot, 'manifest.csv'));

    % --- resume: skip stimuli whose output is already complete ------------
    % Replaces the hardcoded "start at 77", which no longer means anything
    % once naming is decoupled from dir() position.
    nTheta = numel(THETAS);
    todo   = false(numel(stimPaths), 1);
    for i = 1:numel(stimPaths)
        d = fullfile(outRoot, "stim" + stimIDs(i));
        todo(i) = ~(isfolder(d) && numel(dir(fullfile(d, '*.png'))) >= nTheta);
    end

    pathsTodo = stimPaths(todo);
    idsTodo   = stimIDs(todo);
    nTodo     = numel(pathsTodo);

    fprintf('%d stimuli found, %d still to render.\n', numel(stimPaths), nTodo);
    if nTodo == 0; return; end

    % --- create output dirs SERIALLY, before the parfor -------------------
    % mkdir from several workers at once is a race.
    for n = 1:nTodo
        d = fullfile(outRoot, "stim" + idsTodo(n));
        if ~isfolder(d); mkdir(d); end
    end

    if isempty(gcp('nocreate')); parpool; end

    % --- progress reporting (fprintf from a worker goes nowhere useful) ---
    q = parallel.pool.DataQueue;
    nDone = 0;
    afterEach(q, @reportProgress);

    parfor n = 1:nTodo
        processOne(pathsTodo(n), idsTodo(n), outRoot, THETAS, ...
                   APPLY_OUTPUT_GAMMA, USE_COLOR_LUT);
        send(q, idsTodo(n));
    end

    function reportProgress(id)
        nDone = nDone + 1;
        fprintf('%s done (%d/%d)\n', id, nDone, nTodo);
    end
end


% =======================================================================
%  One stimulus: load, composite, convert once, then sweep angles
% =======================================================================
function processOne(fname, stimID, outRoot, thetas, applyGamma, useLUT)

    [rgbIn, ~, alpha] = imread(fname);

    if isempty(alpha)
        error('RotateExplicitPar:noAlpha', '%s has no alpha channel.', fname);
    end

    a    = double(alpha) / 255;
    base = double(rgbIn) / 255;
    base = base .* a + 0.5 .* (1 - a);      % 0.5 = neutral grey background

    [H, W, ~] = size(base);

    % ---- convert to Lab once ---------------------------------------------
    labList = rgb2lab_explicit(reshape(base, [], 3));    % (H*W) x 3

    % ---- collapse to distinct colours ------------------------------------
    if useLUT
        [cols, ~, idx] = unique(labList, 'rows');        % cols: N x 3
    else
        cols = labList;
        idx  = (1:H*W)';
    end

    outDir = fullfile(outRoot, "stim" + stimID);

    L0 = cols(:,1);
    A0 = cols(:,2);
    B0 = cols(:,3);

    for th = thetas
        c = cosd(th);
        s = sind(th);

        rot = [L0, c*A0 - s*B0, s*A0 + c*B0];            % L* untouched

        rgbCols = lab2rgb_explicit(rot, applyGamma);     % N x 3
        out     = reshape(rgbCols(idx,:), H, W, 3);

        fout = fullfile(outDir, sprintf('stim_%s_deg_%03d.png', stimID, th));
        imwrite(out, fout, 'Alpha', a);
    end
end


% =======================================================================
%  STEP 1-3:  sRGB -> Lab   (N x 3 in, N x 3 out)
% =======================================================================
function lab = rgb2lab_explicit(rgb)

    WhitePoint = [0.950456, 1, 1.088754];   % D65, 2 degree observer

    R = invgammacorrection(rgb(:,1));
    G = invgammacorrection(rgb(:,2));
    B = invgammacorrection(rgb(:,3));

    M  = [ 3.240479, -1.537150, -0.498535 ;
          -0.969256,  1.875992,  0.041556 ;
           0.055648, -0.204043,  1.057311 ];
    Mi = inv(M);

    X = Mi(1,1)*R + Mi(1,2)*G + Mi(1,3)*B;
    Y = Mi(2,1)*R + Mi(2,2)*G + Mi(2,3)*B;
    Z = Mi(3,1)*R + Mi(3,2)*G + Mi(3,3)*B;

    fX = labf(X / WhitePoint(1));
    fY = labf(Y / WhitePoint(2));
    fZ = labf(Z / WhitePoint(3));

    lab = [116*fY - 16, 500*(fX - fY), 200*(fY - fZ)];
end


% =======================================================================
%  STEP 5-7:  Lab -> sRGB   (N x 3 in, N x 3 out)
% =======================================================================
function rgb = lab2rgb_explicit(lab, applyGamma)

    WhitePoint = [0.950456, 1, 1.088754];

    L = lab(:,1);
    A = lab(:,2);
    B = lab(:,3);

    fY = (L + 16) / 116;
    fX = fY + A / 500;
    fZ = fY - B / 200;

    X = WhitePoint(1) * labinvf(fX);
    Y = WhitePoint(2) * labinvf(fY);
    Z = WhitePoint(3) * labinvf(fZ);

    M = [ 3.240479, -1.537150, -0.498535 ;
         -0.969256,  1.875992,  0.041556 ;
          0.055648, -0.204043,  1.057311 ];

    R  = M(1,1)*X + M(1,2)*Y + M(1,3)*Z;
    G  = M(2,1)*X + M(2,2)*Y + M(2,3)*Z;
    Bc = M(3,1)*X + M(3,2)*Y + M(3,3)*Z;

    % out-of-gamut: desaturate by adding white, then rescale. Per pixel.
    AddWhite = -min(min(min(R, G), Bc), 0);
    Scale    =  max(max(max(R, G), Bc) + AddWhite, 1);

    R  = (R  + AddWhite) ./ Scale;
    G  = (G  + AddWhite) ./ Scale;
    Bc = (Bc + AddWhite) ./ Scale;

    if applyGamma
        R  = gammacorrection(R);
        G  = gammacorrection(G);
        Bc = gammacorrection(Bc);
    end

    rgb = min(max([R, G, Bc], 0), 1);
end


% =======================================================================
%  Nonlinearities (unchanged, elementwise, shape-agnostic)
% =======================================================================
function fY = labf(Y)
    fY = real(Y .^ (1/3));
    i  = (Y < 0.008856);
    fY(i) = Y(i) * (841/108) + (4/29);
end

function Y = labinvf(fY)
    Y = fY .^ 3;
    i = (Y < 0.008856);
    Y(i) = (fY(i) - 4/29) * (108/841);
end

function R = invgammacorrection(Rp)
    R = real(((Rp + 0.099) / 1.099) .^ (1/0.45));
    i = (R < 0.018);
    R(i) = Rp(i) / 4.5138;
end

function Rp = gammacorrection(R)
    Rp = real(1.099 * R .^ 0.45 - 0.099);
    i  = (R < 0.018);
    Rp(i) = 4.5138 * R(i);
end
