% AUTOMATIC RIB UNFOLDING (UCP-style) — MATLAB reference implementation
% Inspired by: Kolopp et al., "Automatic rib unfolding in postmortem CT" (OpenRib, Unfolded Cylindrical Projection).  % :contentReference[oaicite:1]{index=1}
% Author: <your name/date>
% -------------------------------------------------------------------------

%% USER CONFIG
%dicomDir         = 'PATH/TO/YOUR/DICOM_FOLDER';  % <-- set this
clc; clear;
huBoneThresh     = 150;     % HU threshold to segment cortical/trabecular bone
huBodyThresh     = -300;    % HU threshold to form a torso/body mask
thetaBins        = 720;     % angular sampling (more = smoother tracks)
maxRadiusFactor  = 0.65;    % fraction of body width used as max r for rays
minSlice         = [];      % [] = auto; else restrict (e.g., 50)
maxSlice         = [];      % [] = auto; else restrict (e.g., 300)
savePrefix       = 'rib_unfold';

%% 1) Load DICOM series
%fprintf('Loading DICOM series from: %s\n', dicomDir);
dicom_directory = uigetdir();
fprintf('Loading DICOM series from: %s\n', dicom_directory);
files = dir(fullfile(dicom_directory, '0*'));
if isempty(files); error('No DICOM files in the folder.'); end


% Read metadata to sort slices
info = dicominfo(fullfile(files(1).folder, files(1).name)); % cancle out '.','..'
seriesUID = info.SeriesInstanceUID;
meta = arrayfun(@(f) dicominfo(fullfile(f.folder,f.name)), files(i));
meta = meta(strcmp({meta.SeriesInstanceUID}, seriesUID));
[~, order] = sort(cellfun(@(x) x(3), {meta.ImagePositionPatient}'));
meta = meta(order);
fn = {files.name};
fn = fn(order);

% Read images in order
I = [];
for i = 1:numel(meta)
    I(:,:,i) = double(dicomread(fullfile(files(i).folder, fn{i})));
end

%% . Convert to HU if RescaleSlope/Intercept present
if isfield(meta(1), 'RescaleSlope')
    slope = meta(1).RescaleSlope; intercept = meta(1).RescaleIntercept;
    I = I .* slope + intercept;
end

%% .. Spacing
try
    px = meta(1).PixelSpacing(1);
    py = meta(1).PixelSpacing(2);
catch
    px = 1; py = 1;
end
try
    zLoc  = cellfun(@(m) m.ImagePositionPatient(3), meta);
    dz    = median(abs(diff(zLoc)));
catch
    try
        dz = meta(1).SliceThickness;
    catch
        dz = 1;
    end
end
spacing = [px, py, dz];
fprintf('Volume size: %dx%dx%d  | spacing: [%.3f %.3f %.3f] mm\n', size(I), spacing);

%% 2) Restrict to thorax (optional): pick a z-range heuristically by body area
areaPerSlice = squeeze(sum(I > huBodyThresh, [1 2]));
if isempty(minSlice), minSlice = max(1, find(areaPerSlice > 0.1*max(areaPerSlice), 1, 'first')-5); end
if isempty(maxSlice), maxSlice = min(size(I,3), find(areaPerSlice > 0.1*max(areaPerSlice), 1, 'last')+5); end
I = I(:,:,minSlice:maxSlice);
fprintf('Using slices z = [%d..%d] (count = %d)\n', minSlice, maxSlice, size(I,3));

%% 3) Body & bone masks
fprintf('Segmenting body and bone...\n');
bodyMask = I > huBodyThresh;
% Keep largest 3D component (torso)
bodyMask = bwareaopen(bodyMask, 5000);
CC = bwconncomp(bodyMask, 6);
if CC.NumObjects>1
    vox = cellfun(@numel, CC.PixelIdxList);
    keep = false(CC.NumObjects,1); keep(argmax(vox)) = true;
    bodyMask(:) = false; bodyMask(CC.PixelIdxList{find(keep)}) = true;
end

% Bone mask (cleaned)
boneMask = I > huBoneThresh;
boneMask = boneMask & bodyMask;
for k = 1:size(I,3)
    boneMask(:,:,k) = imclose(boneMask(:,:,k), strel('disk',2));
    boneMask(:,:,k) = bwareaopen(boneMask(:,:,k), 50);
end

%% 4) Per-slice “cage center” via distance transform inside body mask
%    We take the max of the distance map as a robust center (roughly mediastinum/spine region).
fprintf('Estimating per-slice rib-cage centers...\n');
[H, W, Z] = size(I);
centers = zeros(Z,2);
for k = 1:Z
    bm = bodyMask(:,:,k);
    if ~any(bm,'all')
        centers(k,:) = [W/2, H/2];
        continue;
    end
    D = bwdist(~bm);
    [~, idx] = max(D(:));
    [cy, cx] = ind2sub([H, W], idx);
    centers(k,:) = [cx, cy]; % (x,y)
end

% Smooth centers along z for stability
centers = movmean(centers, 7, 1);

%% 5) Build unfolded cylindrical projection (θ × z) with radial MIP over bone mask
fprintf('Building unfolded map (θ × z)...\n');
theta = linspace(-pi, pi, thetaBins+1); theta(end) = []; % [-pi, pi)
% Estimate a max radius per slice from body extent
maxR = zeros(Z,1);
for k = 1:Z
    [yy, xx] = find(bodyMask(:,:,k));
    if isempty(xx)
        maxR(k) = min([W H])/2;
    else
        cx = centers(k,1); cy = centers(k,2);
        maxR(k) = max(sqrt((xx - cx).^2 + (yy - cy).^2));
    end
end
maxR = maxRadiusFactor * maxR;

% Radial sampling step ~ pixel
rStep = 1;

% Prepare outputs
unfoldHU   = nan(thetaBins, Z);    % radial MIP in HU (bones dominate)
unfoldMask = false(thetaBins, Z);  % whether any bone hit on the ray

% For each slice and every angle, sample along a ray and take max HU within bone mask
for k = 1:Z
    cx = centers(k,1); cy = centers(k,2);
    Rk = maxR(k);
    for t = 1:thetaBins
        cth = cos(theta(t)); sth = sin(theta(t));
        rr = 0:rStep:Rk;
        xs = cx + rr .* cth;
        ys = cy + rr .* sth;

        % Bilinear sampling (nearest for masks)
        x0 = round(xs); y0 = round(ys);
        inFOV = x0>=1 & x0<=W & y0>=1 & y0<=H;
        x0 = x0(inFOV); y0 = y0(inFOV);
        idx = sub2ind([H, W], y0, x0);

        if isempty(idx)
            continue;
        end
        boneLine = boneMask(:,:,k);
        hLine    = I(:,:,k);

        maskHits = boneLine(idx);
        if any(maskHits)
            unfoldMask(t,k) = true;
            unfoldHU(t,k)   = max(hLine(idx(maskHits)));
        else
            % fallback: take top HU along ray (gives soft-tissue ribs if threshold misses)
            unfoldHU(t,k) = max(hLine(idx));
        end
    end
end

%% 6) Post-process unfolded map for display
% Normalize to [0,1] for viewing
U = unfoldHU;
% Robust window around bone intensities
pLow  = prctile(U(~isnan(U)), 5);
pHigh = prctile(U(~isnan(U)), 99.5);
U = (U - pLow) / max(1e-6, (pHigh - pLow));
U = min(max(U,0), 1);

% Optional: simple rib-track enhancement (angular 1D top-hat)
Ue = U;
seLen = round(thetaBins/180); % small structuring element along theta
Ue = imtophat(Ue, strel('line', max(3,seLen), 0));

%% 7) Quick visualization
figure('Name','Unfolded ribs (θ × z)','Color','w');
subplot(1,2,1); imshow(U', []); axis on;
title('Unfolded map (normalized HU)'); xlabel('\theta bins (−\pi → \pi)'); ylabel('Axial slice index (z)');
subplot(1,2,2); imshow(Ue', []); axis on;
title('Track-enhanced (top-hat)'); xlabel('\theta bins'); ylabel('z');

% Show a reference axial slice with center
mid = round(Z/2);
figure('Name','Reference axial slice','Color','w');
imshow(windowCT(I(:,:,mid), [pLow pHigh]), []); hold on;
plot(centers(mid,1), centers(mid,2), 'g+', 'MarkerSize', 12, 'LineWidth', 1.5);
title(sprintf('Axial slice %d (center marked)', mid));

%% 8) Optional rib indexing (very simple peak picking per z)
% This is a naive helper to highlight likely rib angles.
fprintf('Estimating rib angle tracks (naive peak picking)...\n');
numRibsToShow = 24; % up to 12 pairs ≈ 24 tracks
ribAngles = nan(numRibsToShow, Z);
for k = 1:Z
    row = Ue(:,k);
    if all(~isfinite(row)), continue; end
    % find peaks in theta
    [pks, locs] = findpeaks(row, 'SortStr','descend', 'NPeaks', numRibsToShow, 'MinPeakProminence', 0.05);
    % sort angularly (left->right)
    locs = sort(locs,'ascend');
    ribAngles(1:min(numel(locs),numRibsToShow), k) = locs(1:min(numel(locs),numRibsToShow));
end

% Overlay tracks
figure('Name','Unfolded with rib-angle overlays','Color','w');
imshow(U', []); hold on;
[rr, zz] = find(~isnan(ribAngles));
plot(ribAngles(~isnan(ribAngles)), zz, 'y.', 'MarkerSize', 4);
title('Unfolded map with naive rib-angle picks');
xlabel('\theta bins'); ylabel('z');

%% 9) Save outputs
outMat = [savePrefix '_unfold.mat'];
save(outMat, 'U', 'Ue', 'unfoldHU', 'unfoldMask', 'centers', 'spacing', 'theta', 'minSlice', 'maxSlice');
fprintf('Saved unfolded data to: %s\n', outMat);

% Helper: simple window function for CT display
function J = windowCT(A, wlwh)
    wl = mean(wlwh); ww = diff(wlwh);
    J = (A - (wl - ww/2)) / ww;
    J = min(max(J,0),1);
end

function i = argmax(v), [~,i] = max(v); end
