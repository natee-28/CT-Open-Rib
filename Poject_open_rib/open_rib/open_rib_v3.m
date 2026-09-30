% AUTOMATIC RIB UNFOLDING (UCP-style) — MATLAB reference implementation
% Inspired by: Kolopp et al., "Automatic rib unfolding in postmortem CT" (OpenRib, Unfolded Cylindrical Projection).  % :contentReference[oaicite:1]{index=1}
% Author: <your name/date>
% -------------------------------------------------------------------------

%% USER CONFIG
%dicomDir         = 'PATH/TO/YOUR/DICOM_FOLDER';  % <-- set this
clc; clear;
huBoneThresh     =  -100;     % HU threshold to segment cortical/trabecular bone
huBodyThresh     =  -150;    % HU threshold to form a torso/body mask
thetaBins        = 1024;     % angular sampling (more = smoother tracks)
maxRadiusFactor  = 1.2;    % fraction of body width used as max r for rays
minSlice         = [];      % [] = auto; else restrict (e.g., 50)
maxSlice         = [];      % [] = auto; else restrict (e.g., 300)
savePrefix       = 'rib_unfold';

%% ..........1: Dicom read Volume......................................
dicom_dir = uigetdir();
%  wait for automatic string of the list 
list = dir(dicom_dir);
% to define which is files not a ".",".."...............................%
if strcmp(list(3).name(1), 'I')
    list_out = sprintf('I*');
elseif strcmp(list(3).name(1), '0')
    list_out = sprintf('0*'); 
end
files = dir(fullfile(dicom_dir,list_out));
% dicom info path 
% info = dicominfo(fullfile(files(1).folder, files(1).name)); % cancle out '.','..'
meta = arrayfun(@(f) dicominfo(fullfile(f.folder,f.name)), files);
%seriesUID = info.SeriesInstanceUID;
seriesUID = meta(1).SeriesInstanceUID;
meta = meta(strcmp({meta.SeriesInstanceUID}, seriesUID));
files = files(strcmp({meta.SeriesInstanceUID},seriesUID));
% ...Compute correct slice order........................................%
iop = meta(1).ImageOrientationPatient;       %row/col directions
row = iop(1:3);
col = iop(4:6);
normal = cross(row,col);         
positions = arrayfun(@(m) dot(m.ImagePositionPatient, normal), meta);
%[~, order] = sort(cellfun(@(x) x(1), {meta.ImagePositionPatient}'));
[~, order] = sort(positions);
meta = meta(order);
files = files(order);
%fn = {files.name};
%fn = fn(order);

%%  Read images in order
I = [];
for i = 1:numel(meta)
    %I(:,:,i) = double(dicomread(fullfile(files(i).folder, fn{i}));
    I(:,:,i) = double(dicomread(fullfile(files(i).folder, files(i).name)));

end
%I = imresize3(I,[512,512,512],'Method','cubic');
%I = permute(fliplr(I),[3,1,2]);


%% . Convert to HU if RescaleSlope/Intercept present
if isfield(meta(1), 'RescaleSlope')
    slope = meta(1).RescaleSlope; intercept = meta(1).RescaleIntercept;
    I = I .* slope + intercept;
end

% Spacing
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
bodyMask = bwareaopen(bodyMask, 100000);
CC = bwconncomp(bodyMask, 6);
if CC.NumObjects>1
    vox = cellfun(@numel, CC.PixelIdxList);
    keep = false(CC.NumObjects,1); keep(argmax(vox)) = true;
    bodyMask(:) = false; bodyMask(CC.PixelIdxList{find(keep)}) = true;
end

% Bone mask (cleaned)
boneMask = I > huBoneThresh;
boneMask = boneMask & bodyMask;
%boneMask = 1.5.*boneMask + bodyMask.*0.5;

for k = 1:size(I,3)
    boneMask(:,:,k) = imclose(boneMask(:,:,k), strel('disk',2));
    boneMask(:,:,k) = bwareaopen(boneMask(:,:,k), 50);
end
boneMask = (bodyMask.*1.5) + (1.0.*boneMask)  ;

%% 4) Per-slice "cage center" via distance transform inside body mask
%    We take the max of the distance map as a robust center (roughly mediastinum/spine region).

fprintf('Estimating per-slice rib-cage centers...\n');
[H, W, Z] = size(I);
centers = zeros(Z,2);
groupSize = round(Z/20);

for gStart = 1:groupSize:Z
    gEnd = min(gStart + groupSize - 1, Z);   % last slice in this group
    k_range = gStart:gEnd;

    % ====== รวมข้อมูลของทั้ง group เข้าเป็น mask เดียว ======
    temp = zeros(H, W);
    for k = k_range
        bm = 1.5 .* bodyMask(:,:,k) + 0.25 .* boneMask(:,:,k);
        temp = temp | bm;   % union ของ mask ทั้ง group
    end
    temp = logical(temp);

    if ~any(temp,'all')
        % ถ้า group นี้ว่าง → ใช้ center กลางภาพ
        cy = H/2;
        cx = W/2;
    else
        % หา largest component ของ group ทั้งก้อน
        CC = bwconncomp(temp);
        stats = regionprops(CC, 'Centroid', 'Area');
        [~, idxMax] = max([stats.Area]);
        cen = stats(idxMax).Centroid;   % [x y]
        cx = cen(1);    % x (columns)
        cy = cen(2);    % y (rows)
    end

    % ====== เซ็ต center เดียวให้ทั้งกลุ่ม ======
    centers(k_range, :) = repmat([cx, cy], numel(k_range), 1);
end


% Smooth centers along z for stability
centers = movmean(centers, 20, 1);
%{
for k=1:Z
centers(k,:) = [W/1.85, H/2.25]; end 
%}
%% 5) Build unfolded cylindrical projection (θ × z) with radial MIP over bone mask

fprintf('Building unfolded map (θ × z)...\n');
theta = linspace((-pi)/2, (3*pi)/2, thetaBins+1); theta(end) = []; % [-pi/2, 3*pi/2)
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
rStep = 1.0;

% Prepare outputs
unfoldHU   = nan(thetaBins, Z);    % radial MIP in HU (bones dominate)
unfoldMask = false(thetaBins, Z);  % whether any bone hit on the ray
U_bone = nan(thetaBins, Z);
U_soft = nan(thetaBins, Z);

soft_k = [];  % just placeholder
%F=zeros(720,1);
soft_tissue = I>-200 & I<250;  
std_xn = @(x) std(nonzeros(x(:)));

% rend= t3_d(I.*boneMask);
% For each slice and every angle, sample along a ray and take max HU within bone mask
for k = 1:Z
    cx = centers(k,1); cy = centers(k,2);
    Rk = maxR(k);
    I_k    = I(:,:,k);
    bone_k = boneMask(:,:,k)>0;
    soft_k = soft_tissue(:,:,k)>0;

    for t = 1:thetaBins
        cth = cos(theta(t)); sth = sin(theta(t));
        rr = 0:rStep:Rk;
        xs = cx + rr .* cth;
        ys = cy + rr .* sth;

        % Bilinear sampling (nearest for masks)
        %x0 = round(xs); y0 = round(ys);
        inFOV = xs>=1 & xs<=W & ys>=1 & ys<=H;
        xs = xs(inFOV); ys = ys(inFOV);
        idx = sub2ind([H, W], ys, xs);

        if isempty(idx)
            continue;
        end
        
        %{
        boneLine = logical(boneMask(:,:,k));
        hLine    = I(:,:,k).*boneMask(:,:,k);%.* soft_tisue(:,:,k) ;
        maskHits =  logical(soft_tisue(idx)+boneLine(idx));%boneLine(idx)+
        %}
        % interpolate HU
        valsHU = interp2(I_k, xs, ys, 'linear', NaN);

        % interpolate masks (nearest is correct for masks)
        valsBone = interp2(double(bone_k), xs, ys, 'nearest', 0) > 0;
        valsSoft = interp2(double(soft_k), xs, ys, 'nearest', 0) > -100;


        if any(valsBone)
          U_bone(t,k) = prctile(valsHU(valsBone), 98);   % 95–99 works well
        end

        % ---- SOFT MAP (between ribs) ----
        softOnly = valsSoft & ~valsBone;
        if any(softOnly)
         U_soft(t,k) = mean(valsHU(softOnly));        % stable soft tissue
        end
        %figure(1); imagesc(unfoldHU.');getframe(gcf);
    end
end
%% ...............show both something..........................................
B = mat2gray(U_bone, [200 1500]);
S = mat2gray(U_soft, [-150 250]);

U_mix = 0.5*B - 0.5*S;

figure(); imshow(flipud(B'),[]); colormap gray; axis image tight;
title('Unfolded Bone + Soft Tissue Blend');

%% 6) Post-process unfolded map for display
% Normalize to [0,1] for viewing
U = U_bone;
% Robust window around bone intensities
pLow  = prctile(U(~isnan(U)), 5);
pHigh = prctile(U(~isnan(U)), 99.5);
U = (U - pLow) / max(1e-6, (pHigh - pLow));
U = min(max(U,0), 1);

% Optional: simple rib-track enhancement (angular 1D top-hat)
Ue = U;
seLen = round(thetaBins/180); % small structuring element along theta
Ue = imtophat(Ue, strel('line', max(16,seLen), 0));
c_Ue = adapthisteq(U, 'ClipLimit', 0.01);

%% 7) Quick visualization
figure('Name','Unfolded ribs (θ × z)','Color','w');
subplot(1,2,1); imshow(U', []); axis on;
title('Unfolded map (normalized HU)'); xlabel('\theta bins (−(pi/2) → (3*pi/2))'); ylabel('Axial slice index (z)');
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
figure(),imshow(imresize(fliplr(flipud(imsharpen(U_bone,Radius=1.25,Amount=0.8).')),[1024,2048],"bilinear"),[])
% Helper: simple window function for CT display
function J = windowCT(A, wlwh)
    wl = mean(wlwh); ww = diff(wlwh);
    J = (A - (wl - ww/2)) / ww;
    J = min(max(J,0),1);
end 

function i = argmax(v), [~,i] = max(v); end



