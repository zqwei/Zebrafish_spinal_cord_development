% code of development video for Shaul's presentation


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
y = y/2;

numLMat = max(preLMatIndex);
LMats    = false(numNeuron, numLMat);
LMat     = preLMat(:, preLMatTime == numTime);
factorIndex = preLMatIndex(preLMatTime == numTime);

for nLMat = 1:numLMat
    nTime = find(preLMatIndex == nLMat, 1, 'first');
    LMats(:, nLMat) = preLMat(:, nTime);
end

figure;

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

plot(y, x, 'ow', 'MarkerFaceColor','w', 'linewidth', 0.01, 'MarkerSize', 10);

for nLMat = 1:numLMat
    plot(y(LMats(:, nLMat)), x(LMats(:, nLMat)), 'ow', 'MarkerFaceColor', mColor(nLMat, :), 'linewidth', 0.01, 'MarkerSize', 10);
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
setPrint(8, 20, 'video', 'pdf')

nStep = 3;
nFactorTime = preLMatTime(find(preLMatIndex==2, 1, 'first'));
LMat = LMats(:, 1);
dffValue  = bsxfun(@plus, zscore(dff(LMat, 1:timePoints(nFactorTime+30)), [], 2), -(1:sum(LMat))'*nStep);
hold on
plot(dffValue', 'Color', mColor(2,:), 'linewidth', 0.5);
gridxy(timePoints(nFactorTime), [], 'color', 'w', 'linestyle', '--')
axis off
setPrint(8, 4, 'Trace', 'pdf')