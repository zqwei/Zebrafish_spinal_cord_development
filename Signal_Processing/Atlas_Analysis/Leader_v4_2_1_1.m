%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial factor correlation analysis
%
% Traces of initial pairs
% Trace feature calculated on inter-event interval distribution 
% Feature similarity calculated by kirmogorov-smirnov statistic
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%
function [totalFactor, followerNonActive, tagList, simIEIList] = Leader_v4_2_1_1(nFile)
% track factor identity using similarity between LMat
addpath('../Func')
setDir;
fileName = fileNames{nFile};
load([tempDatDir '/' fileName '.mat'],'dff', 'activeNeuronMat', 'new_x', 'timePoints');
load([tempDatDir, 'Leader_', fileName, '.mat'], 'leaderPairMetrics', 'activeTime');
% load([tempDatNetDir, fileName, '_spectrogram.mat'], 'spectrogramMatAll');
load([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_2.mat'], 'preLMat', 'preLMatIndex', 'preLMatTime');

newPlotDir = [plotDir  'NeuronPattern_' fileName '/'];

if ~exist([plotDir  'NeuronPattern_' fileName], 'dir')
    mkdir(newPlotDir);
end

spacing    = 0.5;
timeStep   = 600;
followerNonActive = 0;
tagList = [];

bins = 1:30*4; %IEI distribution estimated at frequencies from 0.25s to 30s
prom = 0.25; % peak prominance threshold

simIEIList = cell(0, 1);
figure, 
for nFactor = 1:numel(leaderPairMetrics)
    if nFactor==5
        disp('');
    end
%     r = ceil(nFactor/nCol);
%     c = 4 - (r*nCol - nFactor);
    leaderPairInfo = leaderPairMetrics(nFactor);
    leaderPair      = leaderPairInfo.neuronIndex;
    appearTime      = leaderPairInfo.appearTime;
    firstActiveTime = leaderPairInfo.firstActiveTime;
    preActLevel     = leaderPairInfo.preActLevel;
    if numel(leaderPair) > 2 
        continue;
    end
    clf reset

    % v2: use first activation, appearTime and appearTime+len(pre)
%     tStart = timePoints(max(activeTime(leaderPair)*60))+1; %timePoints(firstActiveTime+1);
    tStart = timePoints(firstActiveTime)+1; %timePoints(firstActiveTime+1);
    tMerge = timePoints(appearTime+1)+1;
    tagNonActive = 0;
    if tMerge-tStart < 5*240 % if merging period shorter than 5 min
        tStart = max(tMerge - 20*240, 1);
        tagNonActive = 1;
        followerNonActive = followerNonActive + 1;
    end
    tEnd   = min(tMerge + tMerge-tStart, timePoints(end));
    t      = min(tStart, timePoints(firstActiveTime+1)+1): tEnd;

    subplot(2, 2, [1, 2]);
    hold on
    timeSeq         = min(t):240:max(t);
    pairCorr        = nan(numel(timeSeq), 1);
    for i = 1:numel(timeSeq)
        tRange      = max(timeSeq(i)-timeStep+1, 1):timeSeq(i)+timeStep;
        pairCorr(i)     = corr(dff(leaderPair(1), tRange)', dff(leaderPair(2), tRange)');
    end
    plot(timeSeq/240, pairCorr, 'g');

    for i = 1:numel(leaderPair)
        yLoc = (i-1)*spacing;
        plot(t/240, dff(leaderPair(i), t)/8+ yLoc);
        text(t(1)/240, yLoc, ['N' num2str(i) ', preActLevel = ' num2str(preActLevel(i), '%.2f')]);
    end
    xlabel('time (min)')
    ylim([-0.2, 1]);
    plot([tMerge, tMerge]/240, [-1, 1], 'k');
    plot([tStart, tStart]/240, [-1, 1], 'k--');
    plot([tEnd, tEnd]/240, [-1, 1], 'k--');
    subplot(2, 2, 3);
    [leaderID, distIEI, numEvents] = plotIEICorr(dff(leaderPair, :), tStart, tMerge, tEnd, prom, bins);
    leaderTag.preActLevel = preActLevel>0.6; 
    leaderTag.mergeLength = (tMerge - tStart)/240;
    leaderTag.events = numEvents;
    IEI = false(2, 1); 

    subplot(2, 2, 4);
    hold on
    mColor = lines(2);
    plot(bins/4, distIEI(:, 1), '-.', 'color', mColor(1, :));
    plot(bins/4, distIEI(:, 2), '-.', 'color', mColor(2, :));
    plot(bins/4, distIEI(:, 3), 'color', mColor(1, :));
    plot(bins/4, distIEI(:, 4), 'color', mColor(2, :));

    legend({'N1-pre', 'N2-pre', 'N1-post', 'N2-post'})
    xlabel('IEI (s)');
    ylabel('cdf');
    hold off

%     simIEIList = [simIEIList; corr(distIEI)];
    simIEIList = [simIEIList; distIEI'*distIEI];
%     % control analysis, time adjacency vs. spectral correlation
%     timeWindowControl = 10*240;
%     tControlStart      = 1:240:t(end)-timeWindowControl;
%     subplot(2, 2, 4);
%     hold on
%     for nCell = 1:2
%         simMetric = nan(0, 2); %timeLag, sigSim
%         maxLag = 0;
%         for t1 = 1:numel(tControlStart)-1
%         sig1 = dff(leaderPair(nCell), tControlStart(t1) + (1:timeWindowControl));
%         iei1 = calculateIEI(sig1, 0.2, bins);
%             for t2 = t1+1:numel(tControlStart)        
%                 sig2 = dff(leaderPair(nCell), tControlStart(t2) + (1:timeWindowControl));
%                 iei2 = calculateIEI(sig2, 0.2, bins);
%                 timeLag = (tControlStart(t2)-tControlStart(t1))/240;
%                 sigSim  = corr(iei1', iei2');
%                 simMetric = [simMetric; [timeLag, sigSim]];
%                 maxLag = max(maxLag, timeLag);
%             end
%         end
% %         scatter(simMetric(:, 1), simMetric(:, 2), '.');
%         [mSim, eSim] = grpstats(simMetric(:, 2), simMetric(:, 1), {'mean', 'std'});
%         errorbar(unique(simMetric(:, 1)), mSim, eSim);
%     end
%     xlabel('time delay');
%     ylabel('IEI auto xcorr');
%     plot([timeWindowControl, timeWindowControl]/240, [-1, 1.5], 'k--');
%     xlim([0, maxLag]);
%     
    if isnan(leaderID)
        titleText = 'Follower non active';
        IEI = [NaN; NaN]; 
    elseif leaderID == 0
        titleText = 'Cannot decide';
        IEI = [0; 0];
    else
        titleText = ['Leader N' num2str(leaderID)];
        IEI(leaderID) = 1; 
    end
    leaderTag.IEICorr = IEI;
    tagList = [tagList; leaderTag];
    subplot(2, 2, [1, 2]);
    title(titleText);
%     setPrint(20, 15, [newPlotDir  'InitialLeaderPairs_IEI_' fileName '_Factor' num2str(nFactor) '_ksdist'], 'pdf');

end
totalFactor = numel(leaderPairMetrics);

end
function [leaderID, distIEI, numEvents] = plotIEICorr(dff, tStart, tMerge, tEnd, prom, bins)
% leaderID = 1 or 2: leader is N1 or N2
% leaderID = NaN: both leader (self correlation)
% leaderID = 0: no leader (mutual correlation)

% raw dff signal, use prom=0.2
sig{1} = dff(1, tStart:tMerge); %1
sig{2} = dff(2, tStart:tMerge); %2
sig{3} = dff(1, tMerge+1:tEnd); %3
sig{4} = dff(2, tMerge+1:tEnd); %4

% % z-scored data, use prom=1 (not good)
% useDFF1 = zscore(dff(1, tStart:tEnd));
% useDFF2 = zscore(dff(2, tStart:tEnd));
% sig{1} = useDFF1(1:tMerge-tStart); 
% sig{2} = useDFF2(1:tMerge-tStart); 
% sig{3} = useDFF1(tMerge-tStart+1:end);
% sig{4} = useDFF2(tMerge-tStart+1:end);

distIEI = nan(numel(bins), 4);
numEvents = zeros(4, 1);

for i = 1:4
    [distIEI(:, i), numEvents(i)] = calculateIEI(sig{i}, prom, bins);
end
% C = corr(distIEI);
% C = distIEI'*distIEI;
if any(numEvents==0)
    leaderID = NaN;
    return;
end
C = nan(4, 4);
for i = 1:4
    for j = i:4
        C(i, j) = max(abs(distIEI(:,  i)-distIEI(:, j))); % ks distance of cdf
    end
end

imagesc(C);
set(gca, 'XTick', 1:4, 'XTickLabel', {'N1-Pre', 'N2-Pre', 'N1-Post', 'N2-Post'})
set(gca, 'YTick', 1:4, 'YTickLabel', {'N1-Pre', 'N2-Pre', 'N1-Post', 'N2-Post'})
colorbar
C(isnan(C)) = nanmin(C(:))-1;
if C(1, 3) >= C(2, 3) && C(1, 4) >= C(2, 4)
    leaderID = 2;
elseif C(2, 3) >= C(1, 3) && C(2, 4) >= C(1, 4)
    leaderID = 1;
else
    leaderID = 0;
end
% C(isnan(C)) = nanmin(C(:))-1; % only consider positive similarity
% if C(1, 3)+C(1,4) > C(2, 3)+C(2, 4)
%     leaderID = 1;
% elseif C(1, 3)+C(1,4) < C(2, 3)+C(2, 4)
%     leaderID = 2;
% else
%     leaderID = NaN;
% end
end

function [iei, nEvents] = calculateIEI(sig, prom, bins)
% inter-event-interval distribution of two signals
% Event detection by local peak detection with threshold
[p, t] = findpeaks(sig,'MinPeakProminence',prom);
% iei = histcounts(t(2:end) - t(1:end-1), bins)./(numel(t)-1);
% if numel(t)>1
%     iei = ksdensity(t(2:end) - t(1:end-1), bins, 'bandwidth', 8);
%     nEvents = numel(t(2:end) - t(1:end-1));
% else
%     iei = zeros(size(bins));
%     nEvents = 0;
% end
usebins = [bins, Inf];
if numel(t)>1
    iei = t(2:end) - t(1:end-1);
    iei(iei>bins(end-1)) = [];
    nEvents = numel(iei);
    iei = histcounts(iei, usebins);
    iei = cumsum(iei, 2)/sum(iei);
else
    iei = zeros(size(bins));
    nEvents = 0;
end

% figure, 
% subplot(2, 1, 1)
% hold on
% plot(1:numel(sig), sig);
% scatter(t, p, 'k');
% hold off
% subplot(2, 1, 2)
% plot(bins, iei);
% disp('');
end