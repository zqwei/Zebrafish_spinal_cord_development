%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Drift Estimation Step 2 - Generate cell cooredinates
% Apply drift estimation onto before and after ablation series
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
function Drift_Estimation_s2(nFile)
addpath('../Func');
setDir;

%4D matrix: nLineage x nTimepoints x 3
load([inputFolderList{nFile, 1} '/DriftCorrection/estimatedDrift.mat'], 'fullDrift');
inputFile = [inputFolderList{nFile, 1}  '\SPM00\TM??????\SPM00_TM??????_CM01_CHN00.klb'];
stackDimensions = size(readImage(recoverFilenameFromPattern(inputFile, 0)));

% drift correction for before
points = dlmread([inputFolderList{nFile, 1} '/Warp/points_before.txt'], ' ');
points(:, 3) = points(:, 3)/6 + 1;
reppoints = repmat(points, [1, 1, numel(fullDrift.before.TP)]);
repdrift = repmat(fullDrift.before.drift, [1, 1, size(points, 1)]);
tracks_smoothed = permute(reppoints, [1, 3, 2]) + permute(repdrift, [3, 1, 2]);

for dim = 1:3
    slicedTrack = squeeze(tracks_smoothed(:, :, dim));
    slicedTrack(slicedTrack<1) = 1;
    slicedTrack(slicedTrack>stackDimensions(dim)) = stackDimensions(dim);
    tracks_smoothed(:, :, dim) = slicedTrack;
end
save([inputFolderList{nFile, 1} '/DriftCorrection/tracks.mat'] , 'tracks_smoothed');

% drift correction for after
points = dlmread([inputFolderList{nFile, 1} '/Warp/points_after.txt'], ' ');
points(:, 3) = points(:, 3)/6 + 1;
reppoints = repmat(points, [1, 1, numel(fullDrift.after.TP)]);
repdrift = repmat(fullDrift.after.drift, [1, 1, size(points, 1)]);
tracks_smoothed = permute(reppoints, [1, 3, 2]) + permute(repdrift, [3, 1, 2]);
for dim = 1:3
    slicedTrack = squeeze(tracks_smoothed(:, :, dim));
    slicedTrack(slicedTrack<1) = 1;
    slicedTrack(slicedTrack>stackDimensions(dim)) = stackDimensions(dim);
    tracks_smoothed(:, :, dim) = slicedTrack;
end
save([inputFolderList{nFile, 3} '/DriftCorrection/tracks.mat'] , 'tracks_smoothed');


