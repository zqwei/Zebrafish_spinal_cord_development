% video for merging



nFile = 3;

addpath('../Func');
setDir;    
fileName          = fileNames{nFile}; %#ok<USENS>    
load([tempDatDir, fileName, '.mat'], 'dff', 'sideSplitter', 'side', 'tracks', 'timePoints', 'new_x', 'new_y', 'new_z'); 
load([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_1.mat'], 'new_activeNeuronMat')
load([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_2.mat'], 'preLMat', 'preLMatIndex', 'preLMatTime');


if ~exist('new_x', 'var'); return; end

x                 = new_x;
y                 = new_y;
z                 = new_z;
numTime           = size(new_activeNeuronMat, 2);
numNeuron         = length(side);
halfActTime       = nan(numNeuron, 1);
timeBin           = 15;
activeThres       = 0.65;

[~, neworder]     = sort(x);
neworder = [neworder(side(neworder)==1); neworder(side(neworder)==2)];
tot_length = length(neworder)+1;

mColor = [       0    0.4470    0.7410
            0.9290    0.6940    0.1250
            0.8500    0.3250    0.0980
            0.4940    0.1840    0.5560
            0.4660    0.6740    0.1880
            0.3010    0.7450    0.9330
            0.6350    0.0780    0.1840
            0.4588    0.4392    0.7020
            0.4000    0.4000    0.4000
            0.9059    0.1608    0.5412]; % 2 and 5 are the finally two color

linew = 2.5;
z = z/max(z) * 1.8;
y = y/2;
for period = 137:-1:130
    timeRange = timePoints(period)+1:timePoints(period)+1200;
    radius = 0.4;

    LMat        = preLMat(:, preLMatTime == period);
    factorIndex = preLMatIndex(preLMatTime == period);
    activeTag   = new_activeNeuronMat(:, period);

    % plot calcium traces
    subplot(1, 8, 1:7);
    hold on
    for i = 1:size(LMat, 1)
        if activeTag(i)
            if ~any(isnan(LMat(i, :))) && sum(LMat(i, :))>0 && size(LMat,2) >= 1
                [~, nFactor] = max(LMat(i, :));
                plot(timeRange, zscore(dff(i, timeRange))+(tot_length - find(neworder==i))*4, 'Color', mColor(factorIndex(nFactor), :), 'linewidth', linew);
            end
        end
    end    
end
xlim([timePoints(130) timePoints(137)+1200])
axis off
set(gca, 'color', 'k')

subplot(1, 8, 8)

period = 152;
LMat        = preLMat(:, preLMatTime == period);
factorIndex = preLMatIndex(preLMatTime == period);
activeTag   = new_activeNeuronMat(:, period);
hold on
if size(LMat,2) >= 1
    for nFactor = 1:size(LMat, 2)
        neuronFactor = LMat(:, nFactor)>0;
        dominateSide = ceil(mean(side(neuronFactor)) - 0.5);
        otherSide    = 3 - dominateSide;
        dominateSideNeuron = neuronFactor & side == dominateSide;
        otherSideNeuron    = neuronFactor & side == otherSide;
        CHPoints = smoothedBoundary(y(dominateSideNeuron), x(dominateSideNeuron), radius);
        plot(CHPoints(:,1), CHPoints(:,2), 'color', mColor(factorIndex(nFactor), :), 'linewid', linew);
        if sum(otherSideNeuron)>0
            if sum(otherSideNeuron) > 0
                CHPoints = smoothedBoundary(y(otherSideNeuron), x(otherSideNeuron), radius);
                plot(CHPoints(:,1), CHPoints(:,2), 'color', mColor(factorIndex(nFactor), :), 'linewid', linew);
            end
        end                

    end
end

period = 130;
LMat        = preLMat(:, preLMatTime == period);
factorIndex = preLMatIndex(preLMatTime == period);
activeTag   = new_activeNeuronMat(:, period);
hold on
plot(y, x, 'ow', 'MarkerFaceColor','w', 'linewidth', 0.01, 'MarkerSize', 10);

for nLMat = 1:size(LMat, 2)
    plot(y(LMat(:, nLMat)), x(LMat(:, nLMat)), 'ow', 'MarkerFaceColor', mColor(factorIndex(nLMat), :), 'linewidth', 0.01, 'MarkerSize', 10);
end

for i = 1:ceil(max(x))-1
    plot([-2, 2], [i, i], '--w');
end
plot([0 0], [0 ceil(max(x))-0.1], '--w')
ylim([0 ceil(max(x))+0.3]);
colorX = linspace(-1, 1, length(mColor)+1);
for nColor = 1:length(mColor)
    plot(colorX(nColor), ceil(max(x))+0.1, 's', 'linewidth', linew, 'color', mColor(nColor,:), 'MarkerFaceColor', mColor(nColor,:), 'MarkerSize', 15)
end
text(-1, ceil(max(x))+0.3, 'Factor Index','color', 'w', 'fontsize', 24)
xlim([-1 1]);
hold off 
axis off
set(gca, 'color', 'k')
set(gca,'YDir','reverse');
set(gcf, 'color', 'k')

setPrint(60, 60, 'Neuron_merge', 'pdf')

timeRange = timePoints(130)+1:timePoints(135);
period = 130;
LMat      = preLMat(:, preLMatTime == period);
global_factor = find(LMat(:, 1));
local_factor  = find(LMat(:, 3));
sum_global_factor = length(global_factor);
dffs = dff([global_factor; local_factor], timeRange);
figure;
imagesc(abs(corr(dffs')), [0 1])


timeRange = timePoints(130)+1:timePoints(135);
period = 130;
LMat      = preLMat(:, preLMatTime == period);
global_factor = find(LMat(:, 1));
local_factor  = find(LMat(:, 3));
sum_global_factor = length(global_factor);
dffs = dff([global_factor; local_factor], timeRange);
figure;
imagesc(abs(corr(dffs')), [0 1])
axis off
setPrint(8, 6, 'Neuron', 'pdf')