%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lineage analysis 1.3
%
% plot the conserved correlation among all 3 datasets
% Use Holm–Bonferroni method of p<0.05
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%

function Lineage_v1_3(datasets)
addpath('../Func');
setDir;

pThres = 0.05;

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

pAll = zeros(numel(leaf), numel(leaf), numel(datasets));
hAll = zeros(numel(leaf), numel(leaf), numel(datasets));

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
    pAll(:, :, nFish) = p;
    hAll(:, :, nFish) = h;
    p_all_positive = p_all_positive + (h>0 & p<pThres);
    p_all_negative = p_all_negative + (h<0 & p<pThres);
end

pAll = sort(pAll, 3, 'descend');
p_corrected = true;
h_count = (sum(hAll>0, 3)==numel(datasets)) - (sum(hAll<0, 3)==numel(datasets));
for i = 1:numel(datasets)
    p_corrected = p_corrected & pAll(:, :, i)<pThres/i;
end
figure, imagesc(triu(p_corrected).*h_count);
set(gca, 'XTick', 1:numel(leaf),'XTickLabel', C(leaf));
set(gca, 'XTickLabelRotation', 45);
set(gca, 'YTick', 1:numel(leaf),'YTickLabel', C(leaf));
set(gca, 'XAxisLocation', 'top');
print( [plotDir, 'LineageCorrMat_conserved_pThres' num2str(pThres) '_HB'], '-dpdf', '-r0');

figure

p_all = p_all_positive - p_all_negative;
p_all(abs(p_all)==1) = 0;


imagesc(triu(p_all, 1));
set(gca, 'XTick', 1:numel(leaf),'XTickLabel', C(leaf));
set(gca, 'XTickLabelRotation', 45);
set(gca, 'YTick', 1:numel(leaf),'YTickLabel', C(leaf));
set(gca, 'XAxisLocation', 'top');
ax1 = imgca;
colorbar(ax1);
print( [plotDir, 'LineageCorrMat_conserved_pThres' num2str(pThres)], '-dpdf', '-r0');