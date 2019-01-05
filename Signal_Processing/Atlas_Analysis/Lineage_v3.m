%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lineage analysis 3: Organize lineage tree to birth stats
%
% 1) trajectory of progenitors
% 2) last division: location, time, mother-children triplets
% 3) lineage tree available
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%
function Lineage_v3(nFile)
addpath('../Func');
setDir;
fileDirName       = fileDirNames{nFile};
fileName          = fileNames{nFile};
load([fileDirName '/', 'dev_data.mat'], 'trackingM', 'leafID');
load([fileDirName '/', 'profile.mat'], 'birthtime');

% trackingM = readMamutXML_trackID_tree(inputXML);
nCells = numel(leafID);
trackignM_splitted = splitTracksbyID(trackingM, leafID, 1:nCells);
trackingM_smoothed = smoothTracks(trackignM_splitted, 1);
birthStats = struct([]);

cellTr = nan(nCells, max(trackingM_smoothed(:, 8))+1, 3); %nCells x nTimepoints x xyz from leaf to root
for i = 1:nCells
    lineageID = trackingM(trackingM(:, 1)==leafID(i), 10);
    if ~isempty(lineageID)
        currentTrack = trackingM_smoothed(trackingM_smoothed(:, 10)==i, :);
        cellTr(i, currentTrack(:, 8)+1, :) = currentTrack(:, 3:5);
    end

    birthStats(i).birthtime     = birthtime(i);  
    birthStats(i).leafNode      = leafID(i);
    birthStats(i).lineageID     = lineageID;
    if ~isnan(birthtime(i))
%         disp(['==== analyzing lineage tree cell #' num2str(i) ' ======'])
        birthLocation = currentTrack(currentTrack(:, 8)==birthtime(i), 3:5);
        lineageTree = trackingM(trackingM(:, 10)==lineageID, :);
        siblingStats = getSibling(lineageTree, leafID(i), birthtime(i));
        birthStats(i).birthLocation = birthLocation;
        birthStats(i).lineageTree  = lineageTree;
        birthStats(i).divisionTriplet  = siblingStats.divisionTriplet;
        birthStats(i).siblingTree = siblingStats.siblingTrackingM;
        birthStats(i).siblingLeafList  = siblingStats.siblingLeafList;
    end
end
  
save([tempDatDir, 'DevelopmentStat_' fileName, '.mat'], 'birthStats', 'cellTr', 'trackingM');
end


