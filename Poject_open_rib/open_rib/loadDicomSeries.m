function [V, metas, files] = loadDicomSeries(root)
% Robust DICOM stack loader (single series in one folder)

if nargin==0 || isempty(root), root = uigetdir(); end
if root==0, error('No folder selected.'); end

% 1) List files, drop folders (this removes '.' and '..') and hidden items
S = dir(root);
S = S(~[S.isdir]);                          % keep only files
S = S(~startsWith({S.name}, '.'));          % drop hidden dotfiles
if isempty(S), error('No files in the folder.'); end

% 2) Keep only DICOM files by probing with dicominfo
keep  = false(numel(S),1);
metas = cell(numel(S),1);
for k = 1:numel(S)
    fn = fullfile(S(k).folder, S(k).name);
    try
        metas{k} = dicominfo(fn);
        keep(k)  = true;
    catch
        % not a DICOM -> ignore
    end
end
S = S(keep); metas = metas(keep);
if isempty(S), error('No DICOM files detected in the selected folder.'); end

% 3) If multiple series are present, pick the dominant SeriesInstanceUID
uids = cellfun(@(m) m.SeriesInstanceUID, metas, 'uni', false);
uid  = char(mode(categorical(uids)));
pick = strcmp(uids, uid);
S = S(pick); metas = metas(pick);

% 4) Determine slice order
haveIPP = cellfun(@(m) isfield(m,'ImagePositionPatient'), metas);
if all(haveIPP)
    % sort along slice normal using ImageOrientationPatient
    iop  = metas{1}.ImageOrientationPatient(:);
    kvec = cross(iop(1:3), iop(4:6));                     % slice-normal
    ipp  = cell2mat(cellfun(@(m)m.ImagePositionPatient(:).', metas, 'uni', false));
    pos  = ipp * kvec;                                    % scalar projection
    [~, order] = sort(pos);
elseif all(cellfun(@(m) isfield(m,'InstanceNumber'), metas))
    [~, order] = sort(cellfun(@(m) double(m.InstanceNumber), metas));
else
    [~, order] = sort([S.datenum]);                       % fallback
end
S = S(order); metas = metas(order);

% 5) Read pixel data
R = double(metas{1}.Rows);  C = double(metas{1}.Columns);  N = numel(S);
V = zeros(R, C, N, 'single');
for k = 1:N
    V(:,:,k) = single(dicomread(fullfile(S(k).folder, S(k).name)));
end

% Optional: apply rescale slope/intercept to get HU if present
slope = isfield(metas{1},'RescaleSlope')     * double(metas{1}.RescaleSlope)     + ~isfield(metas{1},'RescaleSlope');
inter = isfield(metas{1},'RescaleIntercept') * double(metas{1}.RescaleIntercept) + 0*(~isfield(metas{1},'RescaleIntercept'));
V = slope*V + inter;

% return also the ordered file list
files = S;
end
