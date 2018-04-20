%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lineage analysis 3.2: Plot lineage trees
%
% Plot trajectory of progenitors in neural plate space
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%
function Lineage_v3_2(nFile)
addpath('../Func');
setDir;
fileName       = fileNames{nFile};
fileDirName    = fileDirNames{nFile};
load([fileDirName '/', 'dev_data.mat'], 'trackingM', 'leafID');
load([tempDatDir, fileName, '.mat'], 'slicedIndex', 'leafOrder', 'new_x', 'new_y', 'side', 'mnx');
load([tempDatDir, 'DevelopmentStat_' fileName, '.mat'], 'AP', 'LR', 'ori');

leafID = leafID(slicedIndex);
leafID = leafID(leafOrder);

[lineage, annotation] = getTopology(trackingM, leafID);

% calculate new order for lineage trees 
% based on rootNode location on neural plate
rootLocation = nan(numel(lineage), 2); %AP, LR
for i = 1:numel(lineage)
    mamutIDs = lineage(i).MaMuT_ID;
    rootNodeID = mamutIDs(end);
    initialPos = trackingM(trackingM(:, 1)==rootNodeID, 3:5);
    rootLocation(i, :) = (initialPos - ori) * [AP', LR'];
end
[~, newOrder] = sort(rootLocation(:, 2)); % sort by L->R
lineage_sorted = lineage(newOrder);
annotation_sorted = annotation;
for i = 1:size(annotation, 1)
    if ~isnan(annotation(i, 1))
        annotation_sorted(i, 1) = find(newOrder==annotation(i, 1));
    end
end
plotLineageTree(lineage_sorted, annotation_sorted, side, mnx);
setPrint(100, 20, [plotDir 'LineagTreePlot_' fileName '_4_3'], 'pdf');
end

