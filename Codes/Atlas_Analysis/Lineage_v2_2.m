%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lineage analysis 2.2: plot pairwise scatter plot of all datasets
%
% birthtime - activeTime
% birthtime - div angle AP
% neural plate AP - neural plate ML
% neural plate ML - pattern time
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%

function Lineage_v2_2(datasets)
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

    % load and calculate all metrics
    load([tempDatDir, fileName, '.mat'], 'new_x', 'new_y', 'new_z', 'mnx');
    load([tempDatDir, 'Leader_' fileName, '.mat']);
    if nFish == 1
        activeTime(activeTime<0.5) = NaN;
    end
    
    if nFish == 1
        referenceBirth = nanmedian(birthtime);
        referenceAct = nanmedian(patternTime);
        
    else
        birthtime = birthtime - nanmedian(birthtime) + referenceBirth;
        activeTime = activeTime - nanmedian(patternTime) + referenceAct;
        patternTime = patternTime - nanmedian(patternTime) + referenceAct;
    end
    metrics = nan(numel(birthtime), nMetrics+1);
    for i = 1:size(metricList, 1)
        metrics(:, i) = eval(metricList{i, 1});
    end
    metrics(:, end) = nFish;
    metricsAll = [metricsAll; metrics];
end

% birthtime - activeTime
figure,
subplot(1, 4, 1);
gscatter(metricsAll(:, 1), metricsAll(:, 8), metricsAll(:, end));

c = zeros(numel(datasets), 1);
p = zeros(numel(datasets), 1);
for i = 1:numel(datasets)
[c(i), p(i)] = corr(metricsAll(metricsAll(:, end)==i, 1), metricsAll(metricsAll(:, end)==i, 8), 'rows', 'pairwise', 'type','spearman');
disp(['birthtime - activeTime, dataset ' num2str(i) ', c=' num2str(c(i)) ', p=' num2str(p(i))]);
end

xlabel('normalized birthtime (hour)');
ylabel('normalized activation time (hour)');
legend off;
box off;
xlim([0 3]);

% birthtime - div angle AP
subplot(1, 4, 2);
gscatter(metricsAll(:, 1), metricsAll(:, 4), metricsAll(:, end));
for i = 1:numel(datasets)
[c(i), p(i)] = corr(metricsAll(metricsAll(:, end)==i, 1), metricsAll(metricsAll(:, end)==i, 8), 'rows', 'pairwise', 'type','spearman');
disp(['birthtime - div angle AP, dataset ' num2str(i) ', c=' num2str(c(i)) ', p=' num2str(p(i))]);
end
xlabel('normalized birthtime (hour)');
ylabel('division angle AP (degree)');
legend off;
box off;
xlim([0 3]);

% birthtime - patternTime
subplot(1, 4, 3);
gscatter(metricsAll(:, 1), metricsAll(:, 9), metricsAll(:, end));
for i = 1:numel(datasets)
[c(i), p(i)] = corr(metricsAll(metricsAll(:, end)==i, 1), metricsAll(metricsAll(:, end)==i, 8), 'rows', 'pairwise', 'type','spearman');
disp(['birthtime - patternTime ' num2str(i) ', c=' num2str(c(i)) ', p=' num2str(p(i))]);
end
xlabel('normalized birthtime (hour)');
ylabel('normalized pattern time (hour)');
legend off;
box off;
xlim([0 3]);
% neural plate ML - activation time
subplot(1, 4, 4);
gscatter(metricsAll(:, 2), metricsAll(:, 9), metricsAll(:, end));
for i = 1:numel(datasets)
[c(i), p(i)] = corr(metricsAll(metricsAll(:, end)==i, 1), metricsAll(metricsAll(:, end)==i, 8), 'rows', 'pairwise', 'type','spearman');
disp(['neural plate ML - activation time ' num2str(i) ', c=' num2str(c(i)) ', p=' num2str(p(i))]);
end
xlabel('neural plate ML');
ylabel('normalized pattern time (hour)');
legend off;
box off;
xlim([0 1100]);

setPrint(48, 9, [plotDir,  'LineageMetricGroupScatterPlot'], 'pdf');
