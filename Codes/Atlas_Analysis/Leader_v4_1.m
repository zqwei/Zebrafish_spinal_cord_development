%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Track factor identity and analyze the composition of initial factors
%
% Plot size evolution of stable factors
% Plot traces of all initial pairs
%
% Store initial pair metrics to Leader_fileName.mat, including:
% neuronIndex
% activeLength (total # of active windows in merging period)
% appearTime
% firstActiveTime (as defined by activeTime of first neuron)
% xRange
% leaderTag (leaders need to be active > 60% of time in merging period)
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%
function Leader_v4_1(nFile)
% track factor identity using similarity between LMat
addpath('../Func')
setDir;
fileName = fileNames{nFile};
load([tempDatDir '/' fileName '.mat'], 'mnx', 'activeNeuronMat', 'dff', 'timePoints', 'timeStep', 'new_x');
load([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_2.mat'], 'preLMat', 'preLMatIndex', 'preLMatTime');
load([tempDatNetDir, fileName, '_spectrogram.mat'], 'spectrogramMatAll');
load([tempDatDir, 'Leader_', fileName, '.mat'], 'activeTime');

[factorComp, factorSizes] = getFactorIndex(preLMat, preLMatIndex, preLMatTime);
figure, imagesc(factorSizes')
xlabel('Time index');
ylabel('Factor index');
setPrint(8, 6, [plotDir,  'factorSizeIdentityEvol_Network_' fileName], 'pdf');
close

spec = squeeze(mean(spectrogramMatAll, 2));


% leaderPerc = zeros(size(factorSizes, 2), 1);
spacing = 8;

figure
axEnd = zeros(size(factorSizes, 2), 1);
axLen = 0;
ax = [];
leaderPairMetrics = [];
for nFactor = 1:size(factorSizes, 2)
    ax = [ax, subplot(size(factorSizes, 2), 1, nFactor)];
    hold on
    appearTime = find(factorSizes(:, nFactor)>0, 1);
    leaderPair = factorComp{appearTime, nFactor};
    % initial pair: plot traces from first activation to 30 min after appear time
    firstActiveTime = min(activeTime(leaderPair)*60);
    tStart = timePoints(firstActiveTime)+1;
    tEnd   = timePoints(min(appearTime+30, numel(timePoints)))+timeStep;
    activeLength = leaderPair;
    for i = 1:numel(leaderPair)
        activeLength(i) = sum(activeNeuronMat(leaderPair(i), (activeTime(leaderPair(i))*60):appearTime-1));
    end
    
    % leader cell identification
    % 1)merging period >=10 min, 2) active for >=60% in merging period
    activePeriodThres = 10;
    activePercThres   = 0.5;
    if appearTime < activePeriodThres || range(new_x(leaderPair))> 1
        leaderTag = leaderPair * NaN;
%         leaderPairMetrics{nFactor, 1} = NaN;
%         leaderPairMetrics{nFactor, 2} = NaN;
    else
        leaderTag = activeLength/(appearTime-firstActiveTime) > activePercThres;
%         leaderPairMetrics{nFactor, 1} = activeLength;
%         leaderPairMetrics{nFactor, 2} = appearTime-firstActiveTime;
    end
%     leaderPairMetrics{nFactor, 3} = range(new_x(leaderPair));
    for i = 1:numel(leaderPair)
        yLoc = (i-1)*spacing;
        plot(tStart:tEnd, zscore(dff(leaderPair(i), tStart:tEnd))+ yLoc);
        scatter(timePoints(round(activeTime(leaderPair(i))*60))+1, yLoc, 'ok');
        text(tEnd, yLoc, ['Active for ' num2str(activeLength(i), '%02d') ' min, ' num2str(leaderTag(i))]);
    end
    box off
    axis off
    ylim([-1, (i+1)*spacing]);
    axEnd(nFactor) = tEnd;
    axLen = max(axLen, tEnd-tStart);
    plot([timePoints(appearTime), timePoints(appearTime)], [-1, (i+1)*spacing], '--k');
    title(['Pair #' num2str(nFactor) ': ' num2str(firstActiveTime-1) ' - ' num2str(appearTime-1) ' min']);
%     leaderPerc(nFactor) = sum(leaderTag);
    leaderPairInfo.neuronIndex     = leaderPair;
    leaderPairInfo.activeLength    = activeLength;
    leaderPairInfo.appearTime      = appearTime;
    leaderPairInfo.firstActiveTime = firstActiveTime;
    leaderPairInfo.leaderTag       = leaderTag;
    leaderPairInfo.xRange          = range(new_x(leaderPair));
    leaderPairInfo.preActLevel     = activeLength/(appearTime-firstActiveTime);
    leaderPairInfo.prePowerLevel   = mean(spec(leaderPair, tStart:timePoints(appearTime)+1), 2);
    
    leaderPairMetrics = [leaderPairMetrics; leaderPairInfo];
end
for nFactor = 1:size(factorSizes, 2)
    set(ax(nFactor), 'xLim', [axEnd(nFactor)-axLen, axEnd(nFactor)+5000]);
end
setPrint(20, 3*size(factorSizes, 2), [plotDir,  'InitialPairsExample_' fileName], 'pdf');
close

save([tempDatDir, 'Leader_', fileName, '.mat'], 'leaderPairMetrics', '-append');

end

function [factorComp, factorSizes] = getFactorIndex(preLMat, preLMatIndex, preLMatTime)
% factorComp: nTime x nFactor matrix containing cell belongings in each factor
% factorSizes: nFactor x nTime matrix of factor size

numFactors = max(preLMatIndex);
numTime = max(preLMatTime);
factorComp = cell(numTime, numFactors);
factorSizes = zeros(numTime, numFactors);
for nFactor = 1:size(preLMat, 2)
    factorComp{preLMatTime(nFactor), preLMatIndex(nFactor)} = find(preLMat(:, nFactor));
    factorSizes(preLMatTime(nFactor), preLMatIndex(nFactor)) = factorSizes(preLMatTime(nFactor),preLMatIndex(nFactor)) + sum(preLMat(:, nFactor));
end
end