%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extractSignal_fast.m
% fast extraction of signal based on round ROI
% 2.5d method: interpolate the volume to isotropic and get ROI from middle image slice
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%
function extractSignal_fast(parameterDatabase, t)
load(parameterDatabase);
currentStack = readImage(recoverFilenameFromPattern(inputFile, timepoints(t)));
profile = zeros(1, nCells);
for j = 1:nCells
    if any(isnan(tracks_smoothed(j, t, :)))
        profile(j) = NaN;
        continue;
    end
    boundingBox = [round(tracks_smoothed(j, t, 2))-radius, round(tracks_smoothed(j, t, 2))+radius; ...
        round(tracks_smoothed(j, t, 1))-radius, round(tracks_smoothed(j, t, 1))+radius; ...
        floor(tracks_smoothed(j, t, 3)-radius/ratio3D), ceil(tracks_smoothed(j, t, 3) + radius/ratio3D)];
    %prevent bounding box from overflowing image dimensions
    boundingBox(1,1) = max(1, boundingBox(1, 1)); boundingBox(1, 2) = min(size(currentStack, 1), boundingBox(1, 2));
    boundingBox(2,1) = max(1, boundingBox(2, 1)); boundingBox(2, 2) = min(size(currentStack, 2), boundingBox(2, 2)); 
    boundingBox(3,1) = max(1, boundingBox(3, 1)); boundingBox(3, 2) = min(size(currentStack, 3), boundingBox(3, 2));
    
    localWindow = currentStack(boundingBox(1, 1):boundingBox(1, 2), boundingBox(2, 1):boundingBox(2, 2), boundingBox(3, 1):boundingBox(3, 2));
    
    localWindow = imageResize(localWindow, [size(localWindow, 1), size(localWindow, 2), size(localWindow, 3)*ratio3D]);
    windowDimensions = size(localWindow);
    mask = false(2*radius+1, 2*radius+1);
    xRange = cmx;
    yRange = cmy;
    [~, z0] = min(abs(linspace(1, boundingBox(3, 2) - boundingBox(3, 1) + 1, size(localWindow, 3)) - (tracks_smoothed(j, t, 3) - boundingBox(3, 1) + 1)));
    
    invalidIndex = xRange<=0 | xRange>windowDimensions(1) | yRange<=0 | yRange>windowDimensions(2);
    xRange(invalidIndex) = [];
    yRange(invalidIndex) = [];
    mask(sub2ind(windowDimensions, xRange, yRange)) = 1;
    profile(j) = sum(sum(uint16(mask) .* localWindow(:, :, z0)))/numel(xRange);
end

profileName = [outputFolder '\profile.TM' num2str(timepoints(t), '%.6d') '.mat'];
save(profileName, 'profile');