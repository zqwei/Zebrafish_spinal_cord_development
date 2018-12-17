%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lineage analysis 3.3: Plot lineage trees
%
% Plot lineage trees with 
% 1) filtering low-active/incomplete traces
% 2) coloring the bars with EV and halfEV time
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%
function Lineage_v3_3(nFile, visOption)
addpath('../Func');
setDir;
fileName       = fileNames{nFile};
fileDirName    = fileDirNames{nFile};
load([fileDirName '/', 'dev_data.mat'], 'trackingM', 'leafID');
load([tempDatDir, fileName, '.mat'], 'slicedIndex', 'leafOrder', 'new_x', 'new_y', 'side', 'mnx', 'activeNeuronMat');
load([tempDatDir, 'DevelopmentStat_' fileName, '.mat'], 'AP', 'LR', 'ori');
load([tempDatDir, 'Leader_' fileName, '.mat'], 'patternTime');
load([tempDatDir, 'EV_' fileName, '.mat'], 'EVMat');
leafID = leafID(slicedIndex);
leafID = leafID(leafOrder);


% filter out tracks by visualization option
if strcmp(visOption, 'ActMat')
    tagVal = sum(activeNeuronMat, 2) > 0.1*size(activeNeuronMat, 2);
elseif  strcmp(visOption , 'EVMat')
    tagVal = sum(activeNeuronMat, 2) > 0.1*size(activeNeuronMat, 2) & nanmax(EVMat, [], 2)>0.05;
elseif strcmp(visOption , 'patternTime') 
    tagVal = ~isnan(patternTime);
end

tracks2Delete = unique(trackingM(ismember(trackingM(:, 1), leafID(~tagVal)), 10));
tracks2Keep   = unique(trackingM(ismember(trackingM(:, 1), leafID(tagVal)), 10));
trackingM(ismember(trackingM(:,  10), tracks2Delete) & ~ismember(trackingM(:, 10), tracks2Keep), :) = [];
leafID = leafID(tagVal);
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
% [~, newOrder] = sort(abs(rootLocation(:, 2))); % sort by M->L
lineage_sorted = lineage(newOrder);
annotation_sorted = annotation;
for i = 1:size(annotation, 1)
    if ~isnan(annotation(i, 1))
        annotation_sorted(i, 1) = find(newOrder==annotation(i, 1));
    end
end

if strcmp(visOption, 'ActMat')
    plotLineageTreeActivity(lineage_sorted, annotation_sorted, side(tagVal), mnx(tagVal), activeNeuronMat(tagVal, :));
    setPrint(100, 20, [plotDir 'LineagTreePlot' visOption '_' fileName], 'pdf');
elseif strcmp(visOption , 'EVMat')
    % % plot EV tree
    EVMat(isnan(EVMat)) = 0;
    plotLineageTreeEV(lineage_sorted, annotation_sorted, side(tagVal), mnx(tagVal), EVMat(tagVal, :));
    setPrint(100, 20, [plotDir 'LineagTreePlot' visOption '_' fileName], 'pdf');
    EVMat = EVMat./repmat(max(EVMat, [], 2), 1, size(EVMat, 2));
    plotLineageTreeEV(lineage_sorted, annotation_sorted, side(tagVal), mnx(tagVal), EVMat(tagVal, :));
    setPrint(100, 20, [plotDir 'LineagTreePlot' visOption '_' fileName '_normalized'], 'pdf');
elseif strcmp(visOption , 'patternTime') 
    % plot patternTime tree
    ptMat = zeros(size(EVMat));
    pt = round(patternTime*60);
    for i = 1:numel(patternTime)
        if ~isnan(pt(i))
            ptMat(i, pt(i):end) = 1;
        end
    end
    plotLineageTreeActivity(lineage_sorted, annotation_sorted, side(tagVal), mnx(tagVal), ptMat(tagVal, :));
    setPrint(100, 20, [plotDir 'LineagTreePlot' visOption '_' fileName], 'pdf');
end
end

