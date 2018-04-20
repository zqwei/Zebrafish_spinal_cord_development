%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lineage analysis 1.3
%
% plot the conserved correlation among all 3 datasets
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%

function Lineage_v1_3(datasets)
addpath('../Func');
setDir;

% variable name, category (0: dev, 1: func, 2: anatomy, 3: cell type) and
% axis label
metricList = {'birthtime', 0,  'birth time';...
              'abs(neuralPlateLoc(:, 2))', 2, 'neural plate ML';...
              'neuralPlateLoc(:, 1)', 2, 'neural plate AP';...
%               'neuralPlateLoc(:, 3)', 2, 'neural plate DV';...
%               'divAngle(:, 1)', 2, 'division angle radial'; ...
              'divAngle(:, 2)', 2, 'division angle AP'; ...
              'birthLoc(:, 1)', 2, 'birth location AP'; ...
              'abs(birthLoc(:, 2))', 2, 'birth location ML'; ...
%               'birthLoc(:, 3)', 2, 'birth location DV'; ...
              'siblingCloneSize', 0, 'sibling clone size'; ...
              'activeTime', 1, 'activation time'; ...
              'patternTime', 1, 'patterned time'; ...
              'preActLevel', 1, 'pre-pattern activity'; ...
              'new_x', 2, 'AP location'; ...
              'abs(new_y)', 2, 'ML location'; ...
              'new_z', 2, 'DV location'; ...
              'mnx', 3 , 'mnx'; ...
              };
C = metricList(:, 3);
leaf = 1:size(metricList, 1); % leaf order: sorted in 0, 1, 2, 3, within each group grouped in the order of hierachical clustering
p_all_positive = zeros(numel(leaf), numel(leaf));
p_all_negative = zeros(numel(leaf), numel(leaf));




for nFish = 1:numel(datasets)
    nFile = datasets(nFish);
    fileName          = fileNames{nFile};
    
    % load and calculate all metrics
    load([tempDatDir, fileName, '.mat'], 'new_x', 'new_y', 'new_z', 'mnx');
    load([tempDatDir, 'Leader_' fileName, '.mat']);
    
   if nFish == 1
        activeTime(activeTime<0.5) = NaN;
    end
    
    metrics = nan(numel(birthtime), size(metricList, 1));
    for i = 1:size(metricList, 1)
        metrics(:, i) = eval(metricList{i, 1});
    end
    
    [h, p] = corr(metrics(:, leaf), 'rows', 'pairwise', 'type', 'spearman');
    p_all_positive = p_all_positive + (h>0 & p<0.05);
    p_all_negative = p_all_negative + (h<0 & p<0.05);
end

figure

p_all = p_all_positive - p_all_negative;
p_all(abs(p_all)==1) = 0;
% p_all = zeros(numel(leaf), numel(leaf));
% p_all(p_all_positive==numel(datasets)) = 1;
% p_all(p_all_negative==numel(datasets)) = -1;



imagesc(triu(p_all, 1));
set(gca, 'XTick', 1:numel(leaf),'XTickLabel', C(leaf));
set(gca, 'XTickLabelRotation', 45);
set(gca, 'YTick', 1:numel(leaf),'YTickLabel', C(leaf));
set(gca, 'XAxisLocation', 'top');
ax1 = imgca;
% colormap(ax1, colorMap);
colorbar(ax1);
print( [plotDir, 'LineageCorrMat_conserved'], '-dpdf', '-r0');