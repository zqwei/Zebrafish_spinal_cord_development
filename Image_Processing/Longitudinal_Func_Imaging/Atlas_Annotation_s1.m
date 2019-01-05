%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Atlas Annotation Step 1: Convert reference points to landmarks
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
function Atlas_Annotation_s1(inputFolder)
addpath('../Func');
load([inputFolder '/profile/profile.mat'], 'profile_all', 'tracks_smoothed');
markerFile = [inputFolder '/signal/Atlas.tif_resampled.marker'];

ref_points = csvread(markerFile, 1, 0);
ref_points = ref_points(:, [2, 1, 3]);
points = squeeze(tracks_smoothed(:, end, :));
points(:, 3) = points(:, 3) * 6;

ra = ref_points(1:size(ref_points, 1)/2, :);
rb = ref_points(size(ref_points, 1)/2+1:end, :);

figure, scatter3(points(:, 1), points(:, 2), points(:, 3));
hold on, scatter3(ra(:, 1), ra(:, 2), ra(:, 3), 'r');
scatter3(rb(:, 1), rb(:, 2), rb(:, 3), 'g');
hold off

[x, y, z, ~, ~, ~] = convert2atlas3D(points, ra, rb);
save([inputFolder '/tracking/ref_points.mat'], 'points', 'ra', 'rb');
save([inputFolder '/profile/profile.mat'], 'x', 'y', 'z', '-append');