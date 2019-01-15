%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3B
%
% summary composition of initial pairs
% 
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
control_datasets = [3, 4, 7, 12, 10, 15, 16];

addpath('../Func')
setDir;


list = []; % 1activityLength1, 2activityLength2, 3ratio of activeLength, 4difference of activeLength, 5xRange, 6numOfLeaders, 7dataset
for i = 1:numel(control_datasets)
    nFile = control_datasets(i);
    fileName = fileNames{nFile};
    load([tempDatDir, 'Leader_', fileName, '.mat'], 'leaderPairMetrics');    
    for nFactor = 1:numel(leaderPairMetrics)
        activeLength = leaderPairMetrics(nFactor).activeLength;
        xRange = leaderPairMetrics(nFactor).xRange;
        leaderTag = leaderPairMetrics(nFactor).leaderTag;
        activeLength = sort(activeLength);
        if (sum(isnan(leaderTag))==0) && xRange<1
            list = [list; [activeLength(1), activeLength(end), activeLength(1)/activeLength(end), activeLength(end)-activeLength(1), xRange, sum(leaderTag), i]];
        end
    end
end

numLeaders = 0:2;
numLeaderStats = zeros(numel(numLeaders), numel(control_datasets));
 for i = 1:numel(control_datasets)
     numLeaderStats(:, i) = hist(list(list(:, 7)==i, 6), numLeaders);
 end
figure, bar(numLeaders, numLeaderStats, 'stacked');
legend('fish 1', 'fish 2', 'fish 3', 'fish 4', 'fish 5', 'fish 6', 'fish 7')
set(gca, 'XTickLabel', {'N-N', 'L-N', 'L-L'});
ylabel('Number of initial communities');
setPrint(16, 12, [plotDir,  'Figure_3B_InitialLeaderPairs_Type_Summary'], 'pdf');


figure, scatter(list(:, 1), list(:, 2), [], list(:, 7), 'filled');
xLim = ceil(max(max(list(:, 1:2)))/10)*10;
hold on, plot([0, xLim], [0, xLim], 'k--');
xlabel('Activation length neuron 1 (min)');
ylabel('Activation length neuron 2 (min)');
setPrint(16, 12, [plotDir,  'Figure_3B_InitialLeaderPairs_ActivationLength_Summary'], 'pdf');


binRatios = 0:0.2:1;
ratioStats = zeros(numel(binRatios), numel(control_datasets));
 for i = 1:numel(control_datasets)
     ratioStats(:, i) = hist(list(list(:, 7)==i, 3), binRatios);
 end
figure, bar(binRatios, ratioStats, 'stacked');
set(gca, 'XTick', binRatios);
legend('fish 1', 'fish 2', 'fish 3', 'fish 4', 'fish 5', 'fish 6', 'fish 7')
ylabel('Number of initial communities');
xlabel('Ratio of activity level');
setPrint(16, 12, [plotDir,  'Figure_3B InitialLeaderPairs ActivityRatio Summary'], 'pdf');

