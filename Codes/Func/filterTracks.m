function trackingMatrix_out = filterTracks(trackingMatrix, index_ignore)
% filter out the tracking matrix just to keep the entire tracks
% if index_ignore is not empty, also throws away the tracks with track ID
% (trackingMatrix(:, 10) is in index_ignore)

track_id_all = unique(trackingMatrix(:, 10));
track_id_firstTP = trackingMatrix(trackingMatrix(:, 8)==max(trackingMatrix(:, 8)), 10);
brokenTracks = track_id_all(~ismember(track_id_all, track_id_firstTP));
index_ignore = [index_ignore; brokenTracks];

trackingMatrix_out = trackingMatrix;
trackingMatrix_out(ismember(trackingMatrix(:, 10), index_ignore), :) = [];

