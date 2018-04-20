%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5.2, variation of Leader_v5
%
% summary composition of initial pairs using preActLevel
% random simulation for each dataset separately
% both experimental data and simulation shown as mean and std across fish
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
function Leader_v5_2(control_datasets, preActThres)
addpath('../Func')
setDir;

fracExp = zeros(numel(control_datasets), 3); % fraction of NN, NL, LL in each dataset
fracRand = zeros(numel(control_datasets), 3);% expacted fraction of NN, NL, LL in simulated dataset
for i = 1:numel(control_datasets)
    nFile = control_datasets(i);
    fileName = fileNames{nFile};
    load([tempDatDir, 'Leader_', fileName, '.mat'], 'leaderPairMetrics');
    preActExp = []; 
    for nFactor = 1:numel(leaderPairMetrics)        
        xRange = leaderPairMetrics(nFactor).xRange;
        appearTime = leaderPairMetrics(nFactor).appearTime;
        preActLevelPair = leaderPairMetrics(nFactor).preActLevel;
        if (sum(isnan(preActLevelPair))==0) && xRange<1 && numel(preActLevelPair)==2 && appearTime>10
            preActExp = [preActExp; preActLevelPair(1), preActLevelPair(2)];
        end
    end
    fracExp(i, :) = histcounts(sum(preActExp > preActThres, 2), 0:3)/size(preActExp, 1);
    pL = sum(preActExp(:)>preActThres)/numel(preActExp);
    fracRand(i, 1) = pL ^2;
    fracRand(i, 2) = 2*pL*(1-pL);
    fracRand(i, 3) = (1-pL)^2;
end

figure,hold on,
h1 = errorbar(mean(fracExp), std(fracExp), 'r');
h2 = errorbar(mean(fracRand), std(fracRand), 'k');
legend([h1; h2], {'experiment', 'random pairing'});
set(gca, 'XTick', 1:3, 'XTickLabel', {'N-N', 'N-L', 'L-L'});
ylabel('Fraction of initial pairs');
setPrint(8, 6, [plotDir,  'LeaderPairs_avgFish'], 'pdf');

end