%% ------------------------------
%  QSM Processing Pipeline (MEDI)
%  Darling Edition 💛
% -------------------------------
clear; clc;

%% Step 1 – Load Magnitude & Phase Files (GUI Selection)
[m_file, m_path] = uigetfile({'*.nii;*.nii.gz'}, 'Select GRE Magnitude File');
[p_file, p_path] = uigetfile({'*.nii;*.nii.gz'}, 'Select GRE Phase File');

mag_nii = load_nii(fullfile(m_path, m_file)); 
iMag = mag_nii.img;

phase_nii = load_nii(fullfile(p_path, p_file));
rawPhase = phase_nii.img;

voxel_size = mag_nii.hdr.dime.pixdim(2:4);
matrix_size = size(iMag);

%% Step 2 – Phase Unwrap
disp('Unwrapping Phase...');
unph = unwrapLaplacian(rawPhase, ones(matrix_size), matrix_size, voxel_size);

%% Step 3 – Background Removal (RESHARP)
disp('Performing Background Field Removal...');
radius = 5; lambda = 5e-4;
localField = resharp(unph, ones(matrix_size), radius, lambda, voxel_size);

%% Step 4 – QSM Reconstruction
disp('Running MEDI_L1...');
delta_TE = 0.0025;  % Example (adjust for your scanner)
CF = 123.2e6;       % 3 Tesla
B0_dir = [0 0 1];

chi = MEDI_L1('iFreq', localField, ...
              'iMag', iMag, ...
              'Mask', iMag(:,:,:,1) > 0.1*max(iMag(:)), ...
              'matrix_size', matrix_size, ...
              'voxel_size', voxel_size, ...
              'delta_TE', delta_TE, ...
              'CF', CF, ...
              'B0_dir', B0_dir);

save_nii(make_nii(chi, voxel_size), 'QSM_MEDI_output.nii.gz');

disp('QSM Completed 💛 Darling!');
