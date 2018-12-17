%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2
%
% Fit for multi-fish statistics
% 1. number of communities
% 2a. Fraction of non-factored neurons
% 2b. percentage of total active neurons
% 3a. radius of communities
% 3b. size of communities
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
% 

function Figure_2(h, nFile, tag, c)
stats = Figure_2_c(nFile, tag);
totPlots = 6;

figure(h);


% criteria 1: algin by peak factor number
%     xOffset = alignPerc(stats{1}.t, stats{1}.y, 80);
% criteria 2: algin by half recruit time of FracActNeuron
xOffset = stats{2}.fracAct50;

subplot(1, totPlots, 1)
Figure_2_d_1(stats{1}.t-xOffset, stats{1}.y, c)

subplot(1, totPlots, 2)
Figure_2_d_2a(stats{2}.t1-xOffset, stats{2}.y1, c)
Figure_2_d_2b(stats{2}.t2-xOffset, stats{2}.y2, c)

subplot(1, totPlots, 3)
Figure_2_d_3a(stats{3}.t-xOffset, stats{3}.y, c)

subplot(1, totPlots, 4)
Figure_2_d_3b(stats{4}.t-xOffset, stats{4}.y, c)

subplot(1, totPlots, 5)
Figure_2_d_4(stats{5}.t-mean(stats{5}.t)+5, stats{5}.y-xOffset, c, stats{5}.pval<0.05);
view([90, -90]);
set(gca, 'XDir', 'reverse');
subplot(1, totPlots, 6)

Figure_2_d_5(stats{6}.t, stats{6}.y-xOffset, c, stats{6}.pval<0.05);
view([90, -90]);
set(gca, 'XDir', 'reverse');

%     % shift this plot in x so that all points intersect at x=5, y=0
%     x0 = ;
%     t0 = stats{5}.y-xOffset;
%     k = (t0(end)-t0(1))/(x0(end)-x0(1));
%     sOffset = -t0(1)/k+x0(1);
%     Figure_2_d_4(x0+5-sOffset, t0);

addpath('../Func');
setDir;
fileName = fileNames{nFile};
timeOffset = xOffset;
save([tempDatDir, 'Leader_', fileName, '.mat'], 'timeOffset', '-append');

end

%% 1. number of communities
function Figure_2_d_1(t, y, c)
hold on
plot(t, y, 'Color', c);
ylabel('Num factor')
xlabel('Time from peak (hour)')
box off
set(gca, 'TickDir', 'out');
xlim([-2, 5]);
ylim([-1, 7]);
end

%% 2.
% a. Fraction of non-factored neurons
% b. percentage of total active neurons
function Figure_2_d_2a(t, y, c)
hold on
plot(t, y, 'Color', c);
box off
ylabel('Frac active neuron')
xlabel('Normalized time (hour)')
ylim([0, 1.01])
xlim([-2, 5]);
set(gca, 'TickDir', 'out');
end

function Figure_2_d_2b(t, y, c)
hold on
plot(t, y, '-.', 'Color', c);
box off
ylabel('Frac single active neurons')
xlabel('Normalized time (hour)')
ylim([0, 1.01])
xlim([-2, 5]);
set(gca, 'TickDir', 'out');
end

%% 3a. radius of communities
function Figure_2_d_3a(t, y, c)
hold on
plot(t, y, 'Color', c);
box off
xlim([-2, 5]);
ylabel('Radius factor')
xlabel('Normalized time (hour)')
set(gca, 'TickDir', 'out');
end

%% 3b. size of communities
function Figure_2_d_3b(t, y, c)
hold on
plot(t, y, 'Color', c);
ylabel('Size factor')
xlabel('Normalized time (hour)')
box off
ylim([0, 1]);
xlim([-3, 5]);
set(gca, 'TickDir', 'out');
end

%% 4. actTime vs location
function Figure_2_d_4(t, y, c, sig)
hold on
if ~sig
    plot(t, y, '--', 'Color', c);
else
    plot(t, y, 'Color', c);
end
box off
ylabel('Fitted normalized activatoin time (hour)')
xlabel('segment')
set(gca, 'TickDir', 'out');
end

%% 4. mnx vs. median actTime
function Figure_2_d_5(t, y, c, sig)
hold on
if ~sig
    plot(t, y, 'o--', 'Color', c);
else
    plot(t, y, 'o-', 'Color', c);
end
box off
ylabel('Normalized median activation time (hour)')
set(gca, 'TickDir', 'out');
set(gca, 'XTick', 1:2, 'XTickLabel', {'mnx+', 'mnx-'});
xlim([0, 3]);
end

%% related functions
function xOffset = alignPerc(x, y, percVal)
yp     = prctile(y, percVal);
yInd   = find(y >= yp, 1, 'first');
if isempty(yInd); keyboard(); end
xOffset = x(yInd);
end