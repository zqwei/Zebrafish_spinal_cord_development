%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lineage analysis 1: plot correlation matrix, dark background
%
% all metrics prestored in "listLeaderMetrics"
% atlas and mnx information
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%

function Lineage_v1_1(nFile)
addpath('../Func');
setDir;
fileName          = fileNames{nFile};

% load and calculate all metrics
load([tempDatDir, fileName, '.mat'], 'new_x', 'new_y', 'new_z', 'mnx');
load([tempDatDir, 'Leader_' fileName, '.mat']);
diffTime = patternTime - activeTime;

if ~exist('birthtime', 'var')
    return;
end

% variable name and category (0: dev, 1: func, 2: anatomy, 3: cell type)
metricList = {'birthtime', 0; ...
              'activeTime', 1; ...
              'patternTime', 1; ...
              'diffTime', 1; ...
              'factorSize', 1; ...
              'new_x', 2; ...
              'new_y', 2; ...
              'new_z', 2; ...
              'mnx', 3; ...
              'mnxFunc', 3};
metrics = nan(numel(birthtime), size(metricList, 1));

C = metricList(:, 1);
category = cell2mat(metricList(:, 2));

for i = 1:numel(C)
    metrics(:, i) = eval(C{i});
    C{i} = strrep(C{i}, '_', ' ');
end

leaf = 1; % leaf order: sorted in 0, 1, 2, 3, within each group grouped in the order of hierachical clustering
for c = 1:3
    dist = pdist(metrics(:, category==c)', @(x, y)1-corr(x', y', 'row','pairwise'));
    link = linkage(dist);
    order = find(category==c);
    leaf = [leaf; order(optimalleaforder(link, dist))];
end

gamma = 0.6;

nColors = 256;
startColor = [0.2, 0.5, 1];
endColor = [1, 0, 0.2];
midColor = [0, 0, 0];
colorMap = zeros(nColors, 3);
for r = 1:3
    colorMap(:, r) = [linspace(startColor(r), midColor(r), nColors/2)'.^gamma; linspace(midColor(r), endColor(r), nColors/2)'.^gamma];
end

[h, p] = corr(metrics(:, leaf), 'rows', 'pairwise', 'type', 'spearman');
figure
set(0, 'defaultaxeslayer', 'top')
whitebg('white');
imagesc(tril(h, -1), [-.5, .5]);
% im.AlphaData = tril(ones(size(p))) - tril(p);
set(gca, 'XTick', 1:12,'XTickLabel', C(leaf));
set(gca, 'XTickLabelRotation', 45);
set(gca, 'YTick', 1:12,'YTickLabel', C(leaf));
hold on,
[sig_x, sig_y] = find(p<0.05);
text(sig_x(sig_x<sig_y),sig_y(sig_x<sig_y), '*');
ax1 = imgca;
colormap(ax1, colorMap);
colorbar(ax1);
whitebg('black');
set(gcf, 'Color', 'black');
set(gcf, 'InvertHardCopy', 'off', 'PaperPositionMode', 'Auto');
print( [plotDir, 'LineageCorrMat_', fileName], '-dtiff');

figure('Position', [0 , 0, 1000, 300]),
subplot(1, 2, 1)
gscatter(birthtime*1.5/60, activeTime, new_y<0);
xlabel('Birthtime (h)');
ylabel('Active time (h)');
legend('left', 'right', 'Location', 'northwest');
xlim([0, 6]);
subplot(1, 2, 2), 
gscatter(birthtime*1.5/60, patternTime, new_y<0);
hold off
xlim([0, 6]);
xlabel('Birthtime (h)');
ylabel('Pattern time (h)');
legend('left', 'right', 'Location', 'northwest');
whitebg('black');
set(gcf, 'Color', 'black');
set(gcf, 'InvertHardCopy', 'off', 'PaperPositionMode', 'Auto');
print( [plotDir, 'LineageBirthFunc_', fileName], '-dtiff');