function Plot_Segmental_Model(nFile)
% ignore the first and last segments when plotting
addpath('../Func');
setDir;
load([TempDataDir '/tmp_' dataset{nFile} '.mat']);

metrics = {'halfActTime', 'halfEVTime', 'tEV-tAct', 'mnx level', 'factorSize', 'birth time', 'islet'};
for it = 1:7
    validMetrics = 1;
    switch it
        case 1
            lm = fitlm(x, halfActTime, 'linear');
            me = halfActTime - x * lm.Coefficients.Estimate(2) - lm.Coefficients.Estimate(1);
            minR = floor(min(me));
            maxR = ceil(max(me));
            ticks = linspace(minR, maxR, 5);
            ticks = ticks(2:end-1);
        case 2
            lm = fitlm(x, halfEVTime, 'linear');
            me = halfEVTime - x * lm.Coefficients.Estimate(2) - lm.Coefficients.Estimate(1);
            minR = floor(min(me));
            maxR = ceil(max(me));
            ticks = linspace(minR, maxR, 5);
            ticks = ticks(2:end-1);
        case 3
            me = diffHalfTime;
            minR = -0.2;
            maxR = ceil(max(me));
            ticks = linspace(minR, maxR, 5);
            ticks = ticks(2:end-1);
        case 4
            me = mnx_level;
            minR = floor(min(me));
            maxR = ceil(max(me));
            ticks = linspace(minR, maxR, 5);
            ticks = ticks(2:end-1);
        case 5
            factorSize(~mnx) = NaN;
            me = exp(1-factorSize);
            minR = -0.2;
            maxR = 1;
            ticks = 0:0.2:1;
        case 6
            if exist('birthtime', 'var');
                me = birthtime;
                minR = min(birthtime);
                maxR = max(birthtime);
                ticks = linspace(minR, maxR, 5);
                ticks = ticks(2:end-1);
            else
                validMetrics = 0;
            end
        case 7
            if exist('islet', 'var');
                islet(~mnx) = NaN;
                me = islet;
                minR = -0.2;
                maxR = 1;
                ticks = 0:1;
            else
                validMetrics = 0;
            end
    end
    if ~validMetrics
        continue;
    end
    
    me_plot = me(x>=1 & x<=floor(max(x)));
    x_plot = x(x>=1 & x<=floor(max(x)));
    figure,
    mmpolar((x_plot-floor(x_plot))*2*pi, me_plot, 'ob', 'RLim', [minR, maxR], 'TTickDelta', 30, 'RTickValue', ticks, 'TGridColor', [0.8, 0.8, 0.8], 'RGridColor', [0.8, 0.8, 0.8], 'BorderColor', [0.8, 0.8, 0.8]);
    export_fig([PlotDir '/' metrics{it} '_seg_' dataset{nFile} '.pdf'], '-nocrop');
    close
end