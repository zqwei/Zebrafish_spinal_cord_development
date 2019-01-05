%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lineage analysis 1: plot correlation matrix
%
% all metrics prestored in "listLeaderMetrics"
% atlas and mnx information
% (full matrix containing many metics using Lineage_v3_1_1)
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%

function Lineage_v1_2_1(nFile)
addpath('../Func');
setDir;
fileName          = fileNames{nFile};

% load and calculate all metrics
load([tempDatDir, fileName, '.mat'], 'new_x', 'new_y', 'new_z', 'mnx');
load([tempDatDir, 'Leader_' fileName, '.mat']);


if ~exist('birthtime', 'var')
    return;
end


% variable name, category (0: dev, 1: func, 2: anatomy, 3: cell type) and
% axis label
metricList = {'birthtime', 0,  'birth time';...
              'abs(neuralPlateLoc_centroid(:, 2))', 2, 'neural plate ML';...
              'neuralPlateLoc_centroid(:, 1)', 2, 'neural plate AP';...
              'neuralPlateLoc_centroid(:, 3)', 2, 'neural plate DV';...
              'divAngle_centroid(:, 1)', 2, 'division angle radial'; ...
              'divAngle_centroid(:, 2)', 2, 'division angle AP'; ...
              'birthLoc_centroid(:, 1)', 2, 'birth location AP'; ...
              'abs(birthLoc_centroid(:, 2))', 2, 'birth location ML'; ...
              'birthLoc_centroid(:, 3)', 2, 'birth location DV'; ...
              'siblingCloneSize', 0, 'sibling clone size'; ...
              'activeTime', 1, 'activation time'; ...
              'patternTime', 1, 'patterned time'; ...
              'preActLevel', 1, 'pre-pattern activity'; ...
              'new_x', 2, 'AP location'; ...
              'abs(new_y)', 2, 'ML location'; ...
              'new_z', 2, 'DV location'; ...
              'mnx', 3 , 'mnx'; ...
              };
metrics = nan(numel(birthtime), size(metricList, 1));
for i = 1:size(metricList, 1)
    metrics(:, i) = eval(metricList{i, 1});
end
C = metricList(:, 3);
% category = cell2mat(metricList(:, 2));


leaf = 1:size(metricList, 1); % leaf order: sorted in 0, 1, 2, 3, within each group grouped in the order of hierachical clustering


gamma = 0.6;

nColors = 256;
startColor = [0, 0.2, 0.7];
endColor = [0.7, 0, 0];
midColor = [1, 1, 1];
colorMap = zeros(nColors, 3);
for r = 1:3
    colorMap(:, r) = [linspace(startColor(r), midColor(r), nColors/2)'.^gamma; linspace(midColor(r), endColor(r), nColors/2)'.^gamma];
end

[h, p] = corr(metrics(:, leaf), 'rows', 'pairwise', 'type', 'spearman');
figure
set(0, 'defaultaxeslayer', 'top')
whitebg('white');
% h = h(1:end, 2:end);
% p = p(1:end, 2:end);
imagesc(triu(h, 1), [-.5, .5]);
% im.AlphaData = tril(ones(size(p))) - tril(p);
set(gca, 'XTick', 1:numel(leaf),'XTickLabel', C(leaf));
set(gca, 'XTickLabelRotation', 45);
set(gca, 'YTick', 1:numel(leaf),'YTickLabel', C(leaf));
set(gca, 'XAxisLocation', 'top');
hold on,
[sig_x, sig_y] = find(p<0.05);
text(sig_x(sig_x>=sig_y),sig_y(sig_x>=sig_y), '*');
ax1 = imgca;
colormap(ax1, colorMap);
colorbar(ax1);
box off
print( [plotDir, 'LineageCorrMat_', fileName '_full_centroid'], '-dpdf', '-r0');