function trackingM_smoothed = smoothTracks(trackingM, timeStep)

% smooth the input tracking matrix and interpolate to the step of timeStep
% [id, type, x, y, z, radius, parent_id, time, confidence, skeletonId]
trackId = unique(trackingM(:, 10));

trackingM_smoothed = zeros(0, 10);
offset = 0;
for ii = 1:numel(trackId)
    currentTrack = trackingM(trackingM(:, 10)==trackId(ii), :);
    if size(currentTrack, 1)<2
        continue;
    end
    
    currentTPs = currentTrack(:, 8);
    newTPs = min(currentTPs):timeStep:max(currentTPs);
    newTrack = zeros(numel(newTPs), 10);
    newTrack(:, 3:5) = interp1(currentTPs, currentTrack(:, 3:5), newTPs);
    
    newTrack(:, 3) = smooth(newTrack(:, 3), 11, 'rlowess');
    newTrack(:, 4) = smooth(newTrack(:, 4), 11, 'rlowess');
    newTrack(:, 5) = smooth(newTrack(:, 5), 11, 'rlowess');
    newTrack(:, 6) = 10;
    newTrack(:, 1) = offset : offset+numel(newTPs) -1;
    newTrack(1, 7) = -1;
    newTrack(2:end, 7) = newTrack(1:end-1, 1);
    newTrack(:, 8) = newTPs;
    newTrack(:, 10) = trackId(ii);
    offset = offset + numel(newTPs);
    trackingM_smoothed = [trackingM_smoothed; newTrack];
end
