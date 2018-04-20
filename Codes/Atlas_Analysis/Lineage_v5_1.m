%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lineage analysis 5: sibling statistics
%
% 
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%

function Lineage_v5_1(datasets)
addpath('../Func');
setDir;

metricList = {'birthtime/(1.5*60)', 0,  'birth time';... %1
              'abs(neuralPlateLoc(:, 2))', 2, 'neural plate ML';...%2
              'neuralPlateLoc(:, 1)', 2, 'neural plate AP';...%3
              'divAngle(:, 2)', 2, 'division angle AP'; ...%4
              'birthLoc(:, 1)', 2, 'birth location AP'; ...%5
              'birthLoc(:, 2)', 2, 'birth location LR'; ...%6
              'siblingCloneSize', 0, 'sibling clone size'; ...%7
              'activeTime', 1, 'activation time'; ...%8
              'patternTime', 1, 'patterned time'; ...%9
              'preActLevel', 1, 'pre-pattern activity'; ...%10
              'new_x', 2, 'AP location'; ...%11
              'abs(new_y)', 2, 'ML location'; ...%12
              'new_z', 2, 'DV location'; ...%13
              'mnx', 3 , 'mnx'; ...%14
              };
nMetrics = size(metricList, 1);
metricsAll = [];


for nFish = 1:numel(datasets)
    nFile = datasets(nFish);
    fileName          = fileNames{nFile};
    fileDirName       = fileDirNames{nFile};

    % load and calculate all metrics
    load([fileDirName, '/profile.mat'], 'siblings');
    load([tempDatDir, fileName, '.mat'], 'new_x', 'new_y', 'new_z', 'mnx', 'leafOrder', 'slicedIndex');
    load([tempDatDir, 'Leader_' fileName, '.mat']);
    
    ori_id = 1:numel(slicedIndex);
    ori_id = ori_id(slicedIndex);
    ori_id = ori_id(leafOrder);
    for i = 1:size(siblings, 1)
        siblings(i, :) = find(siblings(i, :))==ori_id;
    end
    
    metrics = nan(size(siblings, 1), 2, nMetrics+1);
    for i = 1:size(metricList, 1)
        metric = eval(metricList{i, 1});
        metrics(:, :, i) = metric(siblings);
    end
    metrics(:, end) = nFish;
    metricsAll = [metricsAll; metrics];
    
end

for i =1:nMetircs
    figure,
    gscatter(metricsAll(:, 1, i), metricsAll(:, 2, i));
    title(metricList{i, 3});
end
        