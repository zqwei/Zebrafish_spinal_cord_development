%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial factor correlation analysis - with control of temporal
% continuation
%
% Traces of initial pairs
% Pattern characterized using power spectrum (pmtm + boosting)
% Similarity characterized using correlation - optional version use chi2 distance
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%
function Leader_v4_2_2(nFile)
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

figure, 
for nFactor = 1:numel(leaderPairMetrics)

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
    [leaderID, f, P, C] = plotSpecCorr(dff(leaderPair, :), tStart, tMerge, tEnd);
    titleText = 'Spectral correlation, ';
    leaderTag.preActLevel = preActLevel>0.6; 
    spectralCorr = false(2, 1); 

    if isnan(leaderID)
        titleText = [titleText, 'Both leader'];
        title('Test - Both leader');
        spectralCorr = [1; 1]; 
    elseif leaderID == 0
        titleText = [titleText, 'No leader'];
        title('Test - No leader');
    else
        titleText = [titleText, 'Leader N' num2str(leaderID)];
        spectralCorr(leaderID) = 1; 
        title(['Test - Leader N' num2str(leaderID)]);
    end
    if tagNonActive
        titleText = [titleText, ', follower non-active'];
    end
    
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
%             for t2 = t1+1:numel(tControlStart)        
%                 sig2 = dff(leaderPair(nCell), tControlStart(t2) + (1:timeWindowControl));
%                 timeLag = (tControlStart(t2)-tControlStart(t1))/240;
%                 sigSim  = specCorr(sig1, sig2);
%                 simMetric = [simMetric; [timeLag, sigSim]];
%                 maxLag = max(maxLag, timeLag);
%             end
%         end
%         [mSim, eSim] = grpstats(simMetric(:, 2), simMetric(:, 1), {'mean', 'std'});
%         errorbar(unique(simMetric(:, 1)), mSim, eSim);
%     end
%     xlabel('time delay');
%     ylabel('self spectral correlation');
%     plot([timeWindowControl, timeWindowControl]/240, [0, 1], 'k--');
%     xlim([0, maxLag]);
    subplot(2, 2, 4);
    hold on
%     plot(f, cumsum(P, 2));
    plot(f, P);
    legend({'N1-Pre', 'N2-Pre', 'N1-Post', 'N2-Post'});
    xlabel('frequency (Hz)');
    ylabel('PMTM (Normalized)');
    
    leaderTag.spectralCorr = spectralCorr;
    leaderTag.spectrum = P;
    leaderTag.distM = C;
    leaderTag.freq  = f;
    tagList = [tagList; leaderTag];
    subplot(2, 2, [1, 2]);
    title(titleText);
%     setPrint(20, 15, [newPlotDir  'InitialLeaderPairs_' fileName '_Factor' num2str(nFactor) '_PSD'], 'pdf');
    save([tempDatDir, 'LeaderPairs_PSD', fileName, '.mat'], 'tagList');

end
totalFactor = numel(leaderPairMetrics);

end
function [leaderID, f, P, C] = plotSpecCorr(dff, tStart, tMerge, tEnd)
% leaderID = 1 or 2: leader is N1 or N2
% leaderID = NaN: both leader (self correlation)
% leaderID = 0: no leader (mutual correlation)
    sig1Pre = zscore(dff(1, tStart:tMerge));
    sig2Pre = zscore(dff(2, tStart:tMerge));
    sig1Post = zscore(dff(1, tMerge+1:tEnd));
    sig2Post = zscore(dff(2, tMerge+1:tEnd));
    w = floor(min(numel(sig1Pre), numel(sig1Post))/2)*2;
    fThres = [0, 1]; % signal with frequency above fThres considered as noise
    nw = 101;
%     [Y1pre, f] = periodogram(sig1Pre, [], w, 4);
%     Y2pre      = periodogram(sig2Pre, [], w, 4);
%     Y1post     = periodogram(sig1Post, [], w, 4);
%     Y2post     = periodogram(sig2Post, [], w, 4);
    [Y1pre, f] = pmtm(sig1Pre, nw, w, 4);
    Y2pre      = pmtm(sig2Pre, nw, w, 4);
    Y1post     = pmtm(sig1Post, nw, w, 4);
    Y2post     = pmtm(sig2Post, nw, w, 4);
    fRange     = f<fThres(2) & f>fThres(1);
    P = [Y1pre.*f, Y2pre.*f, Y1post.*f, Y2post.*f];
    P = P(fRange, :)';
    f = f(fRange);
    
    % option 1: use raw PMTM and pearson correlation
    [C, pval] = corr(P');
    leaderID = 0;
    if C(1, 3) > C(2, 3) && C(1, 4) > C(2, 4)
        leaderID = 1;
    elseif C(2, 3) > C(1, 3) && C(2, 4) > C(1, 4)
        leaderID = 2;
    elseif C(1, 3) > C(2, 3) && C(1, 4) < C(2, 4)
        leaderID = NaN;
    end

%     % option 2: normalization of P and similarity using chi-square distance
%     C = nan(4, 4);
%     P = P./repmat(sum(P, 2), 1, sum(fRange));
%     for i = 1:4
%         for j = i:4
%             a = P(i, :);
%             b = P(j, :);
%             C(i, j) = sum((a-b).^2./(a+b))/2; % chi-sq distance
%         end
%     end
%     leaderID = 0;
%     if C(1, 3) < C(2, 3) && C(1, 4) < C(2, 4)
%         leaderID = 1;
%     elseif C(2, 3) < C(1, 3) && C(2, 4) < C(1, 4)
%         leaderID = 2;
%     elseif C(1, 3) < C(2, 3) && C(1, 4) > C(2, 4)
%         leaderID = NaN;
%     end
    
    imagesc(C);
    set(gca, 'XTick', 1:4, 'XTickLabel', {'N1-Pre', 'N2-Pre', 'N1-Post', 'N2-Post'})
    set(gca, 'YTick', 1:4, 'YTickLabel', {'N1-Pre', 'N2-Pre', 'N1-Post', 'N2-Post'})
    colorbar

end
