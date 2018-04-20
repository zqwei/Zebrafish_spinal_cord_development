%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Leader identification based on factor identity tracking
%
% Plot size evolution of stable factors
% Plot traces of all initial pairs
%
% For each neuron calculate
% firstPatternTime
% totalActTime = total active time windows before firstPattern time
% preActLevel = totalActTime/(firstPatternTime - firstActiveTime)
%               firstActiveTime defined per factor
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%
function Leader_v4_2(nFile)
% track factor identity using similarity between LMat
addpath('../Func')
setDir;
fileName = fileNames{nFile};
load([tempDatDir '/' fileName '.mat'], 'mnx', 'activeNeuronMat', 'new_x', 'mnx');
load([tempDatDir, 'Leader_', fileName, '.mat'], 'leaderPairMetrics', 'activeTime');
load([tempDatDir, fileName, '.mat'], 'timePoints');    
load([tempDatNetDir, fileName, '_spectrogram.mat'], 'spectrogramMatAll');
load([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_2.mat'], 'preLMat', 'preLMatIndex', 'preLMatTime');

preActLevel      = nan(size(activeTime));
prePowerLevel    = nan(size(activeTime));
firstPatternTime = nan(size(activeTime));
spec             = squeeze(mean(spectrogramMatAll, 2));
for nCell = 1:numel(activeTime)
    currentLMatIndex        = find(preLMat(nCell, :), 1);
    if ~isempty(currentLMatIndex)
        firstPatternTime = preLMatTime(currentLMatIndex);
        firstActTime            = leaderPairMetrics(preLMatIndex(currentLMatIndex)).firstActiveTime;
        ownActTime              = max(activeTime(nCell)*60, firstActTime);
        if firstPatternTime >= 10 % exclude factors where initial pattern time is very early
            preActLevel(nCell)    = sum(activeNeuronMat(nCell, ownActTime:firstPatternTime-1))/(firstPatternTime-firstActTime);
            prePowerLevel(nCell)  = mean(spec(nCell, timePoints(firstActTime)+1:timePoints(firstPatternTime-1)+1), 2);
        end
    end
end
save([tempDatDir, 'Leader_', fileName, '.mat'], 'preActLevel', 'prePowerLevel', '-append');

preActLevel(isnan(preActLevel)) = [];
bins = 0:0.05:1;
% [f, xi] = ksdensity(preActLevel, 'bandwidth', 0.1);
figure, bar(bins(1:end-1), histcounts(preActLevel, bins)/numel(preActLevel), 'histc');
hold on,
% plot(xi, f, 'r');
xlabel('activity level before factored');
ylabel('percentage');
setPrint(8, 6, [plotDir,  'PreActLevelDistr_' fileName], 'pdf');
% close;

end