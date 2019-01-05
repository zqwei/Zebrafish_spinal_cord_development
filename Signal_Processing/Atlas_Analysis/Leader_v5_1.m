%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Power spectrum analysis and comparison with preActLevel
% 
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
function Leader_v5_1(datasets, preActThres)
addpath('../Func')
setDir;

preActAllPairs = [];
prePowerAllPairs = [];
preActAll = [];
prePowerAll = [];
for i = 1:numel(datasets)
    nFile = datasets(i);
    fileName = fileNames{nFile};
    load([tempDatDir, 'Leader_', fileName, '.mat'], 'preActLevel', 'prePowerLevel', 'leaderPairMetrics', 'activeTime');    
    for nFactor = 1:numel(leaderPairMetrics)        
        xRange = leaderPairMetrics(nFactor).xRange;
        preActLevelPair = leaderPairMetrics(nFactor).preActLevel;
        prePowerPair    = leaderPairMetrics(nFactor).prePowerLevel;
        appearTime      = leaderPairMetrics(nFactor).appearTime;
        if xRange<1 && numel(preActLevelPair)==2 && sum(isnan(prePowerPair))+sum(isnan(preActLevelPair))==0 && appearTime>10
            preActAllPairs = [preActAllPairs; preActLevelPair];
            prePowerAllPairs = [prePowerAllPairs; prePowerPair];
        end
    end
    preActAll = [preActAll; preActLevel];
    prePowerAll = [prePowerAll; prePowerLevel];
end

% % Monte-Carlo simulation for random pairing - baseline estimation
% thresList = 10 .^ (-2.2:0.01:-1);
% nTrials = 10000;
% rankThres = zeros(numel(thresList), 3);
% prePool = prePowerAllPairs(:);
% powerExp = reshape(prePowerAllPairs, 2, [])';
% count = zeros(nTrials, 3);
% for t = 1:numel(thresList)
%     countExp = histcounts(sum(powerExp>thresList(t), 2), 3);
%     for i = 1:nTrials
%         currTrial = reshape(prePool(randperm(numel(prePool))), [], 2);
%         count(i, :) = histcounts(sum(currTrial>thresList(t), 2), 3);
%     end
%     rankThres(t, :) = sum(bsxfun(@gt, countExp, count))/nTrials *100;
% end
% figure, hold on
% plot(thresList, rankThres);
% legend({'pNN', 'pNL', 'pLL'});
% xlabel('Pre-act level threshold');
% ylabel('Pecentile rank of observation');
% targetFun = rankThres(:, 2) - rankThres(:, 1)- rankThres(:, 3); % pNL-pNN-pLL
% prePowerThres = mean(thresList(targetFun==max(targetFun)));
% plot([prePowerThres, prePowerThres], [0, 100], 'k--');
% title(['threshold=' num2str(prePowerThres)]);

validIndex = ~isnan(prePowerAll) & ~isnan(preActAll);
prePowerAll = prePowerAll(validIndex);
preActAll   = preActAll(validIndex);

figure,
subplot(1, 2, 1)
scatter(preActAllPairs, log10(prePowerAllPairs), 'filled');
xlabel('Pre act level')
ylabel('Pre power level (log scale)')
title('Initial pairs only')
subplot(1, 2, 2) 
scatter(preActAll, log10(prePowerAll), 'filled');
xlabel('Pre act level')
ylabel('Pre power level (log scale)')
title('All cells')
setPrint(16, 6, [plotDir,  'PreActLevel_PrePowerLevel'], 'pdf');


figure, 
boxplot(log10(prePowerAll), preActAll>preActThres);
p = ranksum(prePowerAll(preActAll<=preActThres), prePowerAll(preActAll>preActThres));
set(gca, 'XTickLabel', {'Non-Leader', 'Leader'});
ylabel('Pre power level (log scale)')
title(['L vs. N p-value = ' num2str(p)]);
setPrint(8, 6, [plotDir,  'PrePowerLevel_LeaderType'], 'pdf');