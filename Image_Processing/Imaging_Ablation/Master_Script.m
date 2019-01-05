%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Master Script for processing ablation data
% Set configuration files in setDir.m
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org

setDir; % Configure experiment
%% Step 1: Segmentation
for nFile = 1:size(inputFolderList, 1)
    Cell_Segmentation(nFile); % segment cells 
end
% ** Manual Step **
% 1. Curate segmentation result using MTrackJ
% 2. Annotate the motor nerve root in Atlas.tif, save as Atlas.tif_resampled.marker using Vaa3D

%% Step 2: Register time series (dift correction) and before vs. after ablation (warping)
% ** Python Step **
% Run "register_series.py"

%% Step 3: Extract calcium signal
for nFile = 1:size(inputFolderList, 1)
    disp(['Processing fish #' num2str(nFile)]);
    Drift_Estimation_s1(nFile); % Collect drift estimation result
    Drift_Estimation_s2(nFile); % Apply drift estimation to get cell coordinates
    dbList = Extract_Signsl(nFile);
    Atlas_Reconstruction(nFile);
end