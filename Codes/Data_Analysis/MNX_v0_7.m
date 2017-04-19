%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5.  MNX neurons -- active neurons
% 
%  EV and activation as a function of time
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function MNX_v0_7(nFile)            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>    
    mnx               = [];
    load([tempDatDir, fileName, '.mat'], 'mnx', 'new_x'); 
    load([tempDatDir, 'EvoLoading_' fileName, '_v2.mat'], 'networkMat')
    
    if ~exist('mnx', 'var') || ~exist('new_x', 'var')
        return;
    end
    
       
    
    colorSet      = [0.8 0.8 0.8; 0.2 0.2 1.0; 0.2 1.0 0.2;];
        
    numTime       = size(networkMat, 1); %#ok<*NODEF>
    numFactorTime = nan(numTime, 1);
    
    figure; hold on;
    
%     yyaxis left
    
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
        
        
        for nFactor     = find(FactorIndex)
            for mFactor = find(FactorIndex)
                if nFactor < mFactor
                    delayTime = networkCorr.delayMat(nFactor-1, mFactor-1);
                    delayTime = abs(delayTime);
                    delayCorr = abs(networkCorr.corrMat(nFactor-1, mFactor-1));
                    if isnan(delayCorr); delayCorr = 0; end
                    isContra  = networkCorr.IpsiIndex(nFactor-1, mFactor-1) == 0;
                    isOverLap = isOverlapped(new_x, networkTime(nFactor).neuronIndex, networkTime(mFactor).neuronIndex);
                    if isContra
                        markerColor = colorSet((FactorMNX(nFactor) < 1) + (FactorMNX(mFactor) < 1) + 1, :);
                        markerEdgeColor = markerColor;
                        colorRatio  = 1.0;
                        if ~isOverLap
                            markerEdgeColor = [1.0 0.2 0.2]; 
                        end
                        if ~isnan(delayTime)
%                             plot(nTime/60, 28+rand(), 'o', 'markerFaceColor', markerColor, 'markerEdgeColor', markerEdgeColor);
%                         else
%                             markerSize = sizeFASum(networkTime(nFactor).neuronIndex, networkTime(mFactor).neuronIndex);
%                             if markerSize == numUnitsFactorSumMax
%                                 markerSize = 5;
%                             else
%                                 markerSize = 5;
%                             end
%                             markerSize = delayCorr*36;
                            markerSize = 5;
                            plot(nTime/60, delayTime, 'o', 'markerFaceColor', markerColor, 'markerEdgeColor', markerEdgeColor, 'MarkerSize', markerSize);
                        end
                    end
                end
            end
        end
        
    end
    
%     plot([0 numTime/60], [27.5, 27.5], '--k')

%     hold off
    xlabel('Time (hour)')
    ylabel('Contra FA-FA delay (s)')
    xlim([0 numTime/60])
    ylim([0 26])
    box off

%     yyaxis right
    % fit of numFactor curve
    fitResult    = fit((1:numTime)'/60, numFactorTime-2, 'gauss1');
    b            = fitResult.b1;
    a            = fitResult.a1;
    c            = fitResult.c1;
    cr           = c;
    fitResult    = lsqcurvefit(@(p, x) doubleSizedGauss(p, x), [a, b, c, cr], (1:numTime)'/60, numFactorTime);    
    opt1Dim      = doubleSizedGauss(fitResult,(1:numTime)'/60);    
    plot((1:numTime)'/60, opt1Dim./max(opt1Dim)*24,'k-', 'linewid', 2)
    
    
%     hold off
%     xlabel('Time (hour)')
%     ylabel('Num Factors')
%     xlim([0 numTime/60])
%     ylim([0 ceil(max(opt1Dim))])
    box off
    
    setPrint(8*2, 6*2, [plotDir, 'ContraFADelay_' fileName]);
    
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
    
    y = x;
    y(x<=b) = a*exp(-((x(x<=b)-b)/c).^2) + d;
    y(x>b)  = (a+d-dr)*exp(-((x(x>b)-b)/cr).^2) + dr;
end