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
tStart = timePoints(firstActiveTime)+1;
tStart_dff = tStart-1200;
tMerge = timePoints(appearTime);
% tEnd   = timePoints(appearTime+10)+timeStep;
tEnd   = tMerge + tMerge-tStart + 1;

subplot(2, 1, 1)
hold on
offsetPlot = 4;
for i = 1:2
    activeTag = activeNeuronMat(leaderPair(i), :);
    activeTag(1:activeTime(leaderPair(i))*60) = 0;
    plot((tStart_dff:tEnd)/240, dff(leaderPair(i), tStart_dff:tEnd)+ (i-1)*offsetPlot, 'Color', mColor(i, :));
    plot([tStart, tEnd]/240, [(i-1)*offsetPlot - 2, (i-1)*offsetPlot - 2], '--', 'Color', [.8, .8, .8]);
    
    plot(firstActiveTime-5:appearTime+appearTime-firstActiveTime, activeTag(firstActiveTime-5:appearTime+appearTime-firstActiveTime)+ (i-1)*offsetPlot - 2, 'k');
end

box off
% axis off
plot([timePoints(firstActiveTime), timePoints(firstActiveTime)]/240, [-5,10], '--k');
plot([tMerge, tMerge]/240, [-5,10], '--k');

subplot(2, 1, 2)
hold on
tStart_dff   = timePoints(appearTime+40);
offsetPlot = 1;
for i = 1:2
    plot((tStart_dff:tStart_dff+timeStep)/240, dff(leaderPair(i), tStart_dff:tStart_dff+timeStep)+ (i-1)*offsetPlot,'Color', mColor(i, :));
end

set(gcf, 'InvertHardCopy', 'off', 'PaperPositionMode', 'Auto');   
print( [plotDir,  'Figure3C_InitialPair_Example_' fileName '_Pair' num2str(pairID) '_ZoomedView_symwindow.pdf'], '-dpdf', '-r0');

