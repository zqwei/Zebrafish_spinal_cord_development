%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.2. Plot metrics on polar plot
%
% all metrics prestored in "listLeaderMetrics"
% exception: 
% 1) plot halfActTime and halfEVTime first substracting the trend
% 2) plot factorSize in exp scale, and only for mnx+ neurons
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
% 

function Leader_v1_2(nFile) 
addpath('../Func');
setDir;
fileName          = fileNames{nFile};

% calculate factor size
load([tempDatDir, fileName, '.mat'], 'new_x', 'new_y', 'new_z', 'side', 'timePoints', 'mnx');
load([tempDatDir, 'Leader_' fileName, '.mat']);

if ~exist('new_x', 'var')
    return;
end
% calculate other metrics
diffTime = patternTime - activeTime;

% plot atlas with different metrics
x = new_x;
nCol = 4;
nRow = ceil(numel(listLeaderMetrics)/nCol);
figure,
for it = 1:numel(listLeaderMetrics)
    me = eval(listLeaderMetrics{it});
    switch listLeaderMetrics{it}
        case {'activeTime', 'patternTime'}
            lm = fitlm(x, me, 'linear');
            me = me - x * lm.Coefficients.Estimate(2) - lm.Coefficients.Estimate(1);
            minR = floor(min(me));
            maxR = ceil(max(me));
            ticks = linspace(minR, maxR, 5);
            ticks = ticks(2:end-1);
        case 'factorSize'
            me(~mnx) = NaN;
            me = exp(1-me);
            minR = -0.2;
            maxR = 1;
            ticks = 0:0.2:1;
        otherwise
            minR = min(me)-(max(me)-min(me))*0.2;
            maxR = max(me);
            ticks = linspace(minR, maxR, 5);
            ticks = ticks(2:end-1);
    end
    me_plot = me(x>=1 & x<=floor(max(x)));
    x_plot = x(x>=1 & x<=floor(max(x)));

    subplot(nRow, nCol, it)
    mmpolar((x_plot-floor(x_plot))*2*pi, me_plot, 'ob', 'RLim', [minR, maxR], 'TTickDelta', 30, 'RTickValue', ticks, 'TGridColor', [0.8, 0.8, 0.8], 'RGridColor', [0.8, 0.8, 0.8], 'BorderColor', [0.8, 0.8, 0.8]);
    title(listLeaderMetrics{it});
end
    setPrint(8*nCol, 8*nRow, [plotDir,  'LeaderMetricsPolar_' fileName], 'pdf');
    close;
end
