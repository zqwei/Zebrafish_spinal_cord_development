%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5.
%
% summary composition of initial pairs using preActLevel
% Note: use Figure_3B to calculate instead, this is not accurate 
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
function Leader_v5(control_datasets, preActThres)
addpath('../Func')
setDir;

list = []; 
nMulti = 0; %# of pair>2
nXseg = 0; %# of xRange>1
for i = 1:numel(control_datasets)
    nFile = control_datasets(i);
    fileName = fileNames{nFile};
    load([tempDatDir, 'Leader_', fileName, '.mat'], 'leaderPairMetrics');
    for nFactor = 1:numel(leaderPairMetrics)    
        xRange = leaderPairMetrics(nFactor).xRange;
        appearTime = leaderPairMetrics(nFactor).appearTime;
        preActLevelPair = leaderPairMetrics(nFactor).preActLevel;
        if (sum(isnan(preActLevelPair))==0) && appearTime>10
            if xRange<1 && numel(preActLevelPair)==2
                list = [list; preActLevelPair(1), preActLevelPair(2), i];
            elseif numel(preActLevelPair)~=2
                nMulti = nMulti + 1;
            elseif xRange>=1
                nXseg = nXseg + 1;
            end
        end
        
    end
end

% Monte-Carlo simulation for random pairing - baseline estimation
nTrials = 10000;
preActExp = list(:, 1:2);
prePool = preActExp(:);
count = zeros(nTrials, 3);
for i = 1:nTrials
    currTrial = reshape(prePool(randperm(numel(prePool))), [], 2);
    count(i, :) = histcounts(sum(currTrial>preActThres, 2), 3);
end
count_lo = prctile(count, 5);
count_hi = prctile(count, 95);
count_est = mean(count);

numLeaders = 0:2;
numLeaderStats = zeros(3, numel(control_datasets));
 for i = 1:numel(control_datasets)
     numLeaderStats(:, i) = histcounts(sum(preActExp(list(:, 3)==i, :) > preActThres, 2), 0:3);
 end
figure, hold on
h = bar(numLeaders, numLeaderStats, 'stacked', 'EdgeColor', 'none');
h(1).Parent.Parent.Colormap = lines(7);
% legend('fish 1', 'fish 2', 'fish 3', 'fish 4', 'fish 5', 'fish 6', 'fish 7')
set(gca, 'XTick', numLeaders, 'XTickLabel', {'S-S', 'L-S', 'L-L'}, 'YTick', 0:10:40);
ylabel('Number of initial communities');

% fill([numLeaders, fliplr(numLeaders)], [count_lo, fliplr(count_hi)], [.5, .5, .5], 'FaceAlpha', 0.5, 'EdgeColor', 'None');
% plot(numLeaders, count_est, 'k');
errorbar(numLeaders, count_est, count_est-count_lo, count_hi-count_est, 'k');
setPrint(16, 12, [plotDir,  'InitialLeaderPairs_Type_Summary'], 'pdf');
close
% test independence of dataset
p = testIndep(numLeaderStats');
disp(['chi2 test fish independence p=' num2str(p)]);

% try plotting 2 categories
count = [count(:, 2), count(:, 1)+count(:, 3)];
count_lo = prctile(count, 5);
count_hi = prctile(count, 95);
count_est = mean(count);
numLeaderStats = [numLeaderStats(2, :); numLeaderStats(1, :)+numLeaderStats(3, :)];
figure, hold on
h = bar([0 1], numLeaderStats, 'stacked', 'EdgeColor', 'none');
h(1).Parent.Parent.Colormap = lines(7);
% legend('fish 1', 'fish 2', 'fish 3', 'fish 4', 'fish 5', 'fish 6', 'fish 7')
ylabel('Number of initial communities');
bar([2, 3], count_est, 'k');
errorbar([2, 3], count_est, count_est-count_lo, count_hi-count_est, 'k.');
plot([1.5, 1.5], [0 50], 'k--');
set(gca, 'XTick', 0:3, 'XTickLabel', {'L-S', 'L-L or S-S', 'L-S', 'L-L or S-S'}, 'YTick', 0:10:50);
setPrint(16, 12, [plotDir,  'InitialLeaderPairs_Type_2Categories'], 'pdf');



% binRatios = 0:0.2:1;
% ratioStats = zeros(numel(binRatios), numel(control_datasets));
%  for i = 1:numel(control_datasets)
%      ratioStats(:, i) = hist(list(list(:, 3)==i, 1), binRatios);
%  end
% figure, bar(binRatios, ratioStats, 'stacked');
% set(gca, 'XTick', binRatios);
% legend('fish 1', 'fish 2', 'fish 3', 'fish 4', 'fish 5', 'fish 6', 'fish 7')
% ylabel('Number of initial communities');
% xlabel('Ratio of activity level');
% setPrint(16, 12, [plotDir,  'InitialLeaderPairs ActivityRatio Summary'], 'pdf');
% close
end

function p = testIndep(x)
e = sum(x,2)*sum(x)/sum(x(:));
X2 = (x-e).^2./e;
X2 = sum(X2(:));
df = prod(size(x)-[1 1]);
p = 1-chi2cdf(X2,df);
end
