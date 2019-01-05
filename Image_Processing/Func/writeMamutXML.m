function writeMamutXML(trackingMatrix, outputFilename, imageFilename, stackRes)
% write MaMuT xml from trackingMatrix to file
% input:
% trackingMatrix: N x 10 array, where N is the number of spots.
%                 We follow the SWC convention [id, type, x, y, z, radius, parent_id, time, confidence, skeletonId]
% outputFilename: name of xml to export
% imageFilename:  name of the BigDataViewer xml file
% stackRes:       pixel size of [x, y, z]

% header
doc = com.mathworks.xml.XMLUtils.createDocument('TrackMate');
trackMate = doc.getDocumentElement;
trackMate.setAttribute('version', '2.8.1');
model = doc.createElement('Model');
model.setAttribute('spatialunits', 'pixels');
model.setAttribute('timeunits', 'frames');
trackMate.appendChild(model);

% add spots
disp('Adding spots...');
allSpots = doc.createElement('AllSpots');
allTimepoints = unique(trackingMatrix(:, 8));
for i = 1:numel(allTimepoints)
    disp(['    Processing time point ' num2str(i) '/' num2str(numel(allTimepoints))]);
    spotsInFrame = doc.createElement('SpotsInFrame');
    spotsInFrame.setAttribute('frame', num2str(allTimepoints(i)));
    spots = trackingMatrix(trackingMatrix(:, 8)==allTimepoints(i), :);
    for j = 1:size(spots, 1)
        spot = doc.createElement('Spot');
        spot.setAttribute('ID', num2str(spots(j, 1)));
        spot.setAttribute('VISIBILITY', '1');
        if spots(j, 6) < 0
            spot.setAttribute('RADIUS', 'NaN');
        else
            spot.setAttribute('RADIUS', num2str(spots(j, 6), '%.1f'));
        end
        spot.setAttribute('QUALITY', num2str(spots(j, 9), '%.1f')); % confidence
        spot.setAttribute('POSITION_T', num2str(spots(j, 8), '%.1f'));
        spot.setAttribute('POSITION_X', num2str(spots(j, 3)*stackRes(1), '%.13f'));
        spot.setAttribute('POSITION_Y', num2str(spots(j, 4)*stackRes(2), '%.13f'));
        spot.setAttribute('FRAME', num2str(spots(j, 8)));
        spot.setAttribute('POSITION_Z', num2str(spots(j, 5)*stackRes(3), '%.13f'));
        spotsInFrame.appendChild(spot);
    end
    allSpots.appendChild(spotsInFrame);   
end
model.appendChild(allSpots);

% add tracks
disp('Adding tracks...');
allTracks = doc.createElement('AllTracks');
allTrackID = unique(trackingMatrix(:, 10));
for i = 1:numel(allTrackID)
    disp(['    Processing track ' num2str(i) '/' num2str(numel(allTrackID))]);
    track = doc.createElement('Track');
    track.setAttribute('name' , ['Track_' num2str(allTrackID(i))]);
    track.setAttribute('TRACK_ID', num2str(allTrackID(i)));
    edges = trackingMatrix(trackingMatrix(:, 10)==allTrackID(i), :);
    for j = 1:size(edges, 1)
        if edges(j, 7) > 0
            edge = doc.createElement('Edge');
            edge.setAttribute('SPOT_SOURCE_ID', num2str(edges(j, 1)));
            edge.setAttribute('SPOT_TARGET_ID', num2str(edges(j, 7)));
            track.appendChild(edge);
        end
    end
    % if the edge contains a single point, do not add the spot but give a
    % warning instead
    if size(edges, 1) > 1
        allTracks.appendChild(track);
    else
        disp(['Warning: Track ' num2str(allTrackID(i)) ' contains a single spot, omitted'])
    end
end
model.appendChild(allTracks);

% set tracks
filteredTracks = doc.createElement('FilteredTracks');
for i = 1:numel(allTrackID)
    trackID = doc.createElement('TrackID');
    trackID.setAttribute('TRACK_ID', num2str(allTrackID(i)));
    filteredTracks.appendChild(trackID);
end
model.appendChild(filteredTracks);

% footer - BDV xml setting info
disp('Configuring views...');
settings = doc.createElement('Settings');
imageData = doc.createElement('ImageData');
[pathstr,name,ext] = fileparts(imageFilename);
imageData.setAttribute('filename', [name ext]);
imageData.setAttribute('folder', pathstr);
initialFilter = doc.createElement('InitialSpotFilter');
initialFilter.setAttribute('feature', 'QUALITY');
initialFilter.setAttribute('value', '0.0');
initialFilter.setAttribute('idabove', 'true');
settings.appendChild(initialFilter);
settings.appendChild(doc.createElement('SpotFilterCollection'));
settings.appendChild(doc.createElement('TrackFilterCollection'));
settings.appendChild(imageData);
trackMate.appendChild(settings);

% export to file
disp('Writing to file...');
xmlwrite(outputFilename, doc);
end
