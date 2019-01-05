%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.3 Determine the preActLevel threshold for leader vs. non-leader
% Monte-Carlo simulation
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
function preActThres = Leader_v4_3(control_datasets)

addpath('../Func')
setDir;

preActAll = [];
for i = 1:numel(control_datasets)
    nFile = control_datasets(i);
    fileName = fileNames{nFile};
    load([tempDatDir, 'Leader_', fileName, '.mat'], 'preActLevel', 'leaderPairMetrics');
    for nFactor = 1:numel(leaderPairMetrics)        
        xRange = leaderPairMetrics(nFactor).xRange;
        appearTime = leaderPairMetrics(nFactor).appearTime;
        preActLevelPair = leaderPairMetrics(nFactor).preActLevel;
        if (sum(isnan(preActLevelPair))==0) && xRange<1 && numel(preActLevelPair)==2 && appearTime>10
            preActAll = [preActAll; preActLevelPair(1), preActLevelPair(2)];
        end
    end
end

% Monte-Carlo simulation for random pairing - baseline estimation
thresList = 0:0.01:1;
nTrials = 10000;
rankThres = zeros(numel(thresList), 3);
prePool = preActAll(:);
count = zeros(nTrials, 3);
for t = 1:numel(thresList)
    countExp = histcounts(sum(preActAll>thresList(t), 2), 3);
    for i = 1:nTrials
        currTrial = reshape(prePool(randperm(numel(prePool))), [], 2);
        count(i, :) = histcounts(sum(currTrial>thresList(t), 2), 3);
    end
    rankThres(t, :) = sum(bsxfun(@gt, countExp, count))/nTrials *100;
end
figure, hold on
plot(thresList, rankThres);
legend({'pNN', 'pNL', 'pLL'});
xlabel('Pre-act level threshold');
ylabel('Pecentile rank of observation');
targetFun = rankThres(:, 2) - rankThres(:, 1)- rankThres(:, 3); % pNL-pNN-pLL
preActThres = mean(thresList(targetFun==max(targetFun)));
plot([preActThres, preActThres], [0, 100], 'k--');
title(['threshold=' num2str(preActThres)]);
setPrint(8, 6, [plotDir,  'LeaderComposition_MC'], 'pdf');

