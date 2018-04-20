%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3C plot example initial leader pair, preActLevel definition
%
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
addpath('../Func')
setDir;

nFile = 3;
pairID = 5;


fileName = fileNames{nFile};
load([tempDatDir, fileName, '.mat'], 'dff', 'timeStep', 'timePoints', 'activeNeuronMat');
load([tempDatDir, 'Leader_', fileName, '.mat'], 'leaderPairMetrics', 'activeTime');

leaderPair = leaderPairMetrics(pairID).neuronIndex;
appearTime = leaderPairMetrics(pairID).appearTime;
firstActiveTime = leaderPairMetrics(pairID).firstActiveTime;

mColor = lines(numel(leaderPair));

figure, 
tStart_dff = max(timePoints(firstActiveTime)-2400, 1);
tStart = timePoints(firstActiveTime)+1;
tEnd   = timePoints(appearTime+10)+timeStep;

subplot(2, 1, 1)
hold on
offsetPlot = 4;
for i = 1:2
    activeTag = activeNeuronMat(leaderPair(i), :);
    activeTag(1:activeTime(leaderPair(i))*60) = 0;
    plot((tStart_dff:tEnd)/240, dff(leaderPair(i), tStart_dff:tEnd)+ (i-1)*offsetPlot, 'Color', mColor(i, :));
    plot([tStart_dff, tEnd]/240, [(i-1)*offsetPlot - 2, (i-1)*offsetPlot - 2], '--', 'Color', [.8, .8, .8]);
    plot((tStart_dff:timeStep/5:tEnd+1)/240, activeTag(max(1, firstActiveTime-10):appearTime+15)+ (i-1)*offsetPlot - 2, 'k');
end

box off
% axis off
plot([timePoints(firstActiveTime), timePoints(firstActiveTime)]/240, [-5,10], '--k');
plot([timePoints(appearTime), timePoints(appearTime)]/240, [-5,10], '--k');

subplot(2, 1, 2)
hold on
tStart_dff   = timePoints(appearTime+7);
offsetPlot = 1;
for i = 1:2
    plot((tStart_dff:tStart_dff+timeStep)/240, dff(leaderPair(i), tStart_dff:tStart_dff+timeStep)+ (i-1)*offsetPlot,'Color', mColor(i, :));
end

set(gcf, 'InvertHardCopy', 'off', 'PaperPositionMode', 'Auto');   
print( [plotDir,  'Figure3C_InitialPair_Example_' fileName '_Pair' num2str(pairID) '_ZoomedView.pdf'], '-dpdf', '-r0');

