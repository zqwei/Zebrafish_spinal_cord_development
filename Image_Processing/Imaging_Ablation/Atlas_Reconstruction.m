%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Atlas Reconstruction 
% Calculate AP, LR, ML locations based on annotated landmarks
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
function [x, y, z] = Atlas_Reconstruction(nFile)
addpath('../Func');
setDir;

refTimepoint = floor(inputFolderList{nFile, 2}/chn01_interval)*chn01_interval;
load([inputFolderList{nFile, 1} '/Profile/profile.mat'], 'profile_all', 'tracks_smoothed', 'timepoints');
markerFile = [inputFolderList{nFile, 1} '/Markers/Atlas.tif_resampled.marker'];

ref_points = csvread(markerFile, 1, 0);
ref_points = ref_points(:, [2, 1, 3]);
points = squeeze(tracks_smoothed(:, timepoints==refTimepoint, :));
points(:, 3) = points(:, 3) * 6;

ra = ref_points(1:size(ref_points, 1)/2, :);
rb = ref_points(size(ref_points, 1)/2+1:end, :);

% figure, scatter3(points(:, 1), points(:, 2), points(:, 3));
% hold on, scatter3(ra(:, 1), ra(:, 2), ra(:, 3), 'r');
% scatter3(rb(:, 1), rb(:, 2), rb(:, 3), 'g');
% hold off

[x, y, z, ~, ~, ~] = convert2atlas3D(points, ra, rb);
save([inputFolderList{nFile, 1} '/Markers/ref_points.mat'], 'points', 'ra', 'rb');
save([inputFolderList{nFile, 1} '/Profile/profile.mat'], 'x', 'y', 'z', '-append');
save([inputFolderList{nFile, 3} '/Profile/profile.mat'], 'x', 'y', 'z', '-append');