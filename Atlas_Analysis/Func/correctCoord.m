function coordinates_corrected = correctCoord(coordinates)
% coordinates: N x T x 3 matrix representing coordinates in xyzt
% output: 1) drift correct the offset of the entire system
%         2) rotate in xy direction based on PCA of global coordinates

coordinates_corrected = coordinates;
nTime = size(coordinates, 2);
nPoints = size(coordinates, 1);
centers = squeeze(mean(coordinates, 1));

driftVec = centers - repmat(centers(1, :), nTime, 1);

for i = 1:nPoints
    coordinates_corrected(i, :, :) = squeeze(coordinates(i, :, :)) - driftVec;
end

[~, score, ~] = pca(squeeze(reshape(coordinates_corrected(:, :, 1:2), nTime*nPoints, 2)));
coordinates_corrected(:, :, 1:2) = reshape(score, nPoints, nTime, 2) + 2000;
end