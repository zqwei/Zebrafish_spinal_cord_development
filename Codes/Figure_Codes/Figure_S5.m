%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure S5 activation and pattern time
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%
addpath('../Func')
setDir;
load([fileDirNames{3} '/profile.mat'], 'timepoints');
load([tempDatDir '/' fileNames{3} '.mat'], 'dff', 'activeNeuronMat');
load([tempDatDir '/EV_' fileNames{3} '.mat'], 'EVMat');

nNeuron = 68;
mColor  = lines(2);
numTime           = size(EVMat, 2);

figure,
subplot(3, 1, 1); plot((1:size(dff, 2))/(3600*4), dff(nNeuron, :), 'k');
xlim([0 numTime/60]);
ylim([-0.3, 1.3])
box off
ylabel('dff');

subplot(3, 1, 2:3);
set(gca, 'XTick', 0:4);
timeIndex = ~isnan(EVMat(nNeuron, :));
actCurrNeuron   = activeNeuronMat(nNeuron, :);
smoothActIndex  = smooth(double(actCurrNeuron), timeBin);
if sum(actCurrNeuron(1:9)) < 5
    smoothActIndex(1:5) = 0;
end
firstActiveTime = find(smoothActIndex > activeThres, 1, 'first');
if ~isempty(firstActiveTime)
    removeIndex     = find(actCurrNeuron==0);
    removeIndex(removeIndex<firstActiveTime) = [];
    removeTimeIndex = false(1, numTime);
    removeTimeIndex(removeIndex) = true; % remove all zeros after first activation time
    timeIndex = timeIndex & ~ removeTimeIndex;
else
    removeTimeIndex = false(1, numTime);
    timeIndex       = timeIndex & ~ removeTimeIndex;
end

if ~isempty(firstActiveTime)
    mdl = fitglm(find(~removeTimeIndex)/60, activeNeuronMat(nNeuron, ~removeTimeIndex), 'linear', 'Distribution', 'binomial');
    beta = mdl.Coefficients.Estimate;
    halfActTime = -beta(1)/beta(2);
    if (mean(activeNeuronMat(nNeuron, 1:ceil(numTime*.5))) > 0.80 ||  mean(activeNeuronMat(nNeuron, :)) > 0.80) ...
            && halfActTime<0 || halfActTime>numTime/60
        halfActTime = 0;
    end
    if (halfActTime<0 || halfActTime>numTime/60); halfActTime = nan; end
end

if sum(timeIndex) > 10 && ~isempty(firstActiveTime) % && halfActTime(nNeuron)>=0
    paddingLength            = 20;
    zerosPadding             = zeros(1, paddingLength+1);
    init_params              = [0, quantile(EVMat(nNeuron, timeIndex), 0.95), halfActTime, 1];
    
    if halfActTime > numTime || isnan(halfActTime)
        init_params(3)       = firstActiveTime/60;
        init_params(4)       = 0.001;
    end
    
    [fitParams, fitResult]   = sigm_fit([(-paddingLength:0)/60, find(timeIndex)/60], [zerosPadding, EVMat(nNeuron, timeIndex)], [0, nan, nan, nan], init_params, false);
    ypred                    = fitResult.ypred(paddingLength+2: end);
    ypredlowerCI             = fitResult.ypredlowerCI(paddingLength+2:end);
    ypredupperCI             = fitResult.ypredupperCI(paddingLength+2:end);
    RSquare(nNeuron)         = 1 - mean((EVMat(nNeuron, timeIndex)' - ypred).^2)./var(EVMat(nNeuron, timeIndex)'); %#ok<UDIM>
    halfEVTime     = fitParams(3);
    %             if nNeuron==2
    %                 keyboard();
    %             end
    halfEVTimeCI(nNeuron, :)    = fitResult.paramCI(2, :);
    if halfEVTime < 0 || halfEVTime > numTime/60
        halfEVTime  = nan;
    end
    %             [~, fitResult]           = sigm_fit(find(timeIndex)/60, EVMat(nNeuron, timeIndex)');
    %             RSquare(nNeuron)         = 1 - mean((EVMat(nNeuron, timeIndex)' - fitResult.ypred).^2)./var(EVMat(nNeuron, timeIndex)');
end

[hAx, hl_EV, hl_act] = plotyy(find(timeIndex)/60, ypred, find(~removeTimeIndex)/60, nan(size(mdl.Fitted.Probability)));
hl_EV.Color = mColor(1, :);
hl_act.Color = mColor(2, :);
hAx(2).YAxis.Color = mColor(1, :);
hAx(1).YAxis.Color = mColor(2, :);
% plot EV curve
hold(hAx(1));
plot(hAx(1), find(timeIndex)/60, EVMat(n, timeIndex), 'o', 'Color', mColor(1, :));
plot(hAx(1),find(timeIndex)/60, ypred, '-', 'linewid', 2.0, 'Color', mColor(1, :));
plot(hAx(1),find(timeIndex)/60, ypredlowerCI, '-', 'linewid', 0.5, 'Color', mColor(1, :));
plot(hAx(1), find(timeIndex)/60, ypredupperCI, '-', 'linewid', 0.5, 'Color', mColor(1, :));
plot([halfEVTime, halfEVTime], [0, 0.5],':',  'linewid', 2.0, 'Color', mColor(1, :));
xlim(hAx(1), [0 numTime/60]);
ylim(hAx(1), [0 1]);
ylabel(hAx(1), 'EV');

% plot activation curve
hold(hAx(2));
plot(hAx(2), find(~removeTimeIndex)/60, actCurrNeuron(~removeTimeIndex), 'd', 'Color', mColor(2, :));
% plot(hAx(2), find(~removeTimeIndex)/60, mdl.Fitted.Probability, '-', 'linewid', 2.0, 'Color', mColor(2, :));
% plot(hAx(2), [halfActTime, halfActTime], [0, 0.5],':', 'linewid', 2.0, 'Color', mColor(2, :));
plot(hAx(2), [firstActiveTime/60, firstActiveTime/60], [0, 1],':', 'linewid', 2.0, 'Color', mColor(2, :));
xlim(hAx(2), [0 numTime/60]);
ylim(hAx(2), [0 1]);
ylabel(hAx(2), 'active state');
hold off;
ylim([0 1])
xlim([0 numTime/60])
xlabel('Time (hour)')
box off

setPrint(16, 12, [plotDir,  'Figure_S5 Act_EV_example'], 'pdf');

