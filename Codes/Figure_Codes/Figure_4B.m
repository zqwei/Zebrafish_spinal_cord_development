%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Figure 4B: plot MNX neurons -- LR phase delay
%
% adapted from MNX_v0_7
% Finite delay -- circle
% Circle size  -- factor size
% 180419: add random phase percentage at top,  use hpfilt to fit #Factors
% Yinan Wan
% figure: run Figure_4B(3)
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function Figure_4B(nFile)
addpath('../Func');
setDir;
fileName          = fileNames{nFile}; %#ok<USENS>
load([tempDatDir, 'FALONO_', fileName, '.mat'], 'dumpDuplicatedFactorLONOM');

load([tempDatDir, fileName, '.mat'], 'mnx', 'new_x');
if ~exist([tempDatDir, 'EvoLoading_' fileName, '_v2.mat'], 'file'); return; end
load([tempDatDir, 'EvoLoading_' fileName, '_v2.mat'], 'networkMat')

if ~exist('mnx', 'var') || ~exist('new_x', 'var')
    return;
end

colorSet      = [0.5 0.5 0.5; 0.7461 0.3828 0.6367; 0.3516 0.7305 0.9141];

numTime       = size(networkMat, 1); %#ok<*NODEF>
numFactorTime = nan(numTime, 1);
figure; 
subplot(4, 1, 2:4)
hold on;

% fit of numFactor curve
tmpLONOM = [zeros(100, 1); dumpDuplicatedFactorLONOM; ones(100, 1)*2];
lambda  = l1tf_lambdamax(tmpLONOM);
opt1Dim = hpfilter(tmpLONOM, lambda*2);
opt1Dim  = opt1Dim(101:end-100);
[hAx, h1, h2] =  plotyy((1:numTime)'/60, nan(size(opt1Dim)), (1:numTime)'/60, opt1Dim);
h1.Color = [0 0 0];
h2.Color = [1 0.2 0.2];
hAx(2).YAxis.Color = [1 0.2 0.2];
hAx(1).YAxis.Color = [0 0 0];
set(h2, 'LineWidth', 2);

xlabel('Time (hour)')
ylabel('Num Factors')
xlim([0 numTime/60])
ylim([0 ceil(max(opt1Dim))])
box off

fracRandom = nan(numTime, 1);
for nTime     = 1:numTime
    networkTime    = networkMat{nTime, 1}; %#ok<*USENS>
    networkCorr    = networkMat{nTime, 2};
    numUnitsFactor = [networkTime.nUnit];
    FactorIndex    = numUnitsFactor > 1;
    FactorIndex(1) = false;
    numFactorTime(nTime) = sum(FactorIndex);
    FactorMNX      = [networkTime.mnxFrac];
    if sum(FactorIndex) == 0
        continue;
    end
    
    %         numUnitsFactor = numUnitsFactor(FactorIndex);
    %         numUnitsFactorSum = numUnitsFactor + numUnitsFactor';
    %         numUnitsFactorSum(eye(length(numUnitsFactor))==1) = 0;
    %         numUnitsFactorSumMax = max(numUnitsFactorSum(:));
    
    nRandom  = 0;
    nNonRand = 0;
    for nFactor     = find(FactorIndex)
        for mFactor = find(FactorIndex)
            if nFactor < mFactor
                delayTime = networkCorr.delayMat(nFactor-1, mFactor-1);
                delayTime = abs(delayTime);
                delayCorr = abs(networkCorr.corrMat(nFactor-1, mFactor-1));
                if isnan(delayCorr); delayCorr = 0; end
                isContra  = networkCorr.IpsiIndex(nFactor-1, mFactor-1) == 0;
%                 isOverLap = isOverlapped(new_x, networkTime(nFactor).neuronIndex, networkTime(mFactor).neuronIndex);
                if isContra
                    markerColor = colorSet((FactorMNX(nFactor) < 1) + (FactorMNX(mFactor) < 1) + 1, :);
%                     markerEdgeColor = markerColor;
%                     colorRatio  = 1.0;
%                     if ~isOverLap
%                         markerEdgeColor = [1.0 0.2 0.2]; % clusters directly oposite to each other
%                     end
                    if ~isnan(delayTime) && delayTime<10
                        %                             plot(nTime/60, 28+rand(), 'o', 'markerFaceColor', markerColor, 'markerEdgeColor', markerEdgeColor);
                        %                         else
                        %                             markerSize = sizeFASum(networkTime(nFactor).neuronIndex, networkTime(mFactor).neuronIndex);
                        %                             if markerSize == numUnitsFactorSumMax
                        %                                 markerSize = 5;
                        %                             else
                        %                                 markerSize = 5;
                        %                             end
                        markerSize = (delayCorr > 0.3)*10 + 5;
                        %                             markerSize = delayCorr*36;
                        %                             markerSize = 5;
%                         plot(hAx(1), nTime/60, delayTime, 'o', 'markerFaceColor', markerColor, 'markerEdgeColor', markerEdgeColor, 'MarkerSize', markerSize);
                        plot(hAx(1), nTime/60, delayTime, 'o', 'markerFaceColor', markerColor, 'markerEdgeColor', 'none', 'MarkerSize', markerSize);
                        nNonRand = nNonRand + 1;

                    else
                        %                         markerSize = (delayCorr > 0.3)*10 + 5;
                        nRandom = nRandom + 1;
%                         plot(nTime/60, nRandom + 11, 's', 'markerFaceColor', markerColor, 'markerEdgeColor', markerEdgeColor, 'MarkerSize', markerSize);
                    end
                end
            end
        end
    end
    fracRandom(nTime) = nRandom/(nRandom+nNonRand);
end

%     plot([0 numTime/60], [27.5, 27.5], '--k')

%     hold off
xlabel('Time (hour)')
ylabel(hAx(1),'Contra FA-FA delay (s)')
ylabel(hAx(2),'Number of Factor')
xlim([0 numTime/60])
set(hAx(1), 'YTick', 0:2:10, 'YLim', [0, 12]);
set(hAx(1), 'TickDir', 'out');
set(hAx(2), 'TickDir', 'out');
box off
hold off

subplot(4, 1, 1)
bar((1:numTime)/60, fracRandom*100, 'k');
ylabel('random phase (%)');
set(gca, 'TickDir', 'out');
box off
xlim([0 numTime/60])


setPrint(8*2, 6*2, [plotDir, 'ContraFADelay_' fileName '_randomportion'], 'pdf');

end

function y = sizeFAMax(factor1, factor2)
y = max(length(factor1), length(factor2));
end

function y = sizeFAMin(factor1, factor2)
y = min(length(factor1), length(factor2));
end

function y = sizeFASum(factor1, factor2)
y = sum([length(factor1), length(factor2)]);
end


function y = isOverlapped(xLoc, factor1, factor2)
xLoc1  = xLoc(factor1);
xLoc2  = xLoc(factor2);
min1   = min(xLoc1);
max1   = max(xLoc1);
min2   = min(xLoc2);
max2   = max(xLoc2);

y      = (min1 > max2) | (min2 > max1);
y      = ~y;
end

function y = doubleSizedGauss(p, x)
p      =  num2cell(p);

[a, b, c, cr] = deal(p{:});
d = 0;
dr = 2;
%     dr = 0;

y = x;
y(x<=b) = a*exp(-((x(x<=b)-b)/c).^2) + d;
y(x>b)  = (a+d-dr)*exp(-((x(x>b)-b)/cr).^2) + dr;
end
