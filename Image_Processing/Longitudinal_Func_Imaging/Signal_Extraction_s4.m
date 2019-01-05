%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Signal Extraction Step 4: Interpolate and smooth curated tracks
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%
function Signal_Extraction_s4(inputFolder, ratio3D, rlowessWindow)
addpath('../Func');
inputTracks = [inputFolder '.simpleTrack\signal\data-mamut.xml'];
outputTracks = [inputFolder '.simpleTrack\tracking\tracks_smoothed.mat'];
outputFormat = 1;     % outputFormat 1 - 4D matrix, 2 - swc, 3 - mamutxml
writeIncomplete = 0;  % 0 - write only complete tracks, 1 - write also incomplete tracks
dataXML = [inputFolder '.simpleTrack\signal\data.xml']; % provide only if outputFormat=3

disp('Reading MaMuT xml file...')
swc = readMamutXML_trackID(inputTracks);

% translate MaMuT coordinate to pixel space
swc(:, 5) = swc(:, 5)/ratio3D;
swc(:, [3, 4, 5]) = swc(:, [4, 3, 5]) + 1; 


lineages = unique(swc(:, 10));
timePoints = unique(swc(:, 8));
tracks_smoothed = nan(numel(lineages), numel(timePoints), 3);
disp('interpolating and smoothing node location ...')
for i = 1:numel(lineages)
    currentTrack = swc(swc(:, 10)==lineages(i), :);
    timeRange = ismember(timePoints, currentTrack(:, 8));
    % interpolate the skipped time points
    tracks_smoothed(i, :, :) = interp1(currentTrack(:, 8), currentTrack(:, 3:5), timePoints);
    % rlowess smoothing of the 3 dimensions with the specified span
    for j = 1:3
        tracks_smoothed(i, timeRange, 1) = smooth(currentTrack(:, 3), rlowessWindow, 'rlowess');
        tracks_smoothed(i, timeRange, 2) = smooth(currentTrack(:, 4), rlowessWindow, 'rlowess');
        tracks_smoothed(i, timeRange, 3) = smooth(currentTrack(:, 5), rlowessWindow, 'rlowess');
    end
end

if ~writeIncomplete
    disp('killing incomplete tracks')
    incomplete_flag = false(numel(lineages));
    for i = 1:numel(lineages)
        if (any(isnan(tracks_smoothed(i, :, 1))))
            incomplete_flag(i) = true;
        end
    end
    tracks_smoothed(incomplete_flag, :, :) = [];
    lineages(incomplete_flag) = [];
end

if outputFormat == 1
    save(outputTracks, 'tracks_smoothed', 'lineages', 'timePoints');
else
    swc_updated = export4DMatrix2SWC(tracks_smoothed, lineages, timePoints);
    if outputFormat == 2
        save(outputTracks, 'swc_updated');
    elseif outputFormat == 3
        % map pixel space back to MaMuT coordinates
        swc_updated(:, [4, 3, 5]) = swc_updated(:, [3, 4, 5]) - 1 ;
        writeMamutXML(swc_updated, outputTracks, dataXML, [1 1 ratio3D]);
    end
end