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
    
    figure; hold on;
    
    for nTime     = 1:numTime
        networkTime    = networkMat{nTime, 1}; %#ok<*USENS>
        networkCorr    = networkMat{nTime, 2};
        numUnitsFactor = [networkTime.nUnit];
        FactorIndex    = numUnitsFactor > 1;
        FactorIndex(1) = false;
        FactorMNX      = [networkTime.mnxFrac];
        if sum(FactorIndex) == 0
            continue;
        end
        
        for nFactor     = find(FactorIndex)
            for mFactor = find(FactorIndex)
                if nFactor < mFactor
                    delayTime = networkCorr.delayMat(nFactor-1, mFactor-1);
                    delayTime = abs(delayTime);
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
                            plot(nTime/60, delayTime, 'o', 'markerFaceColor', markerColor, 'markerEdgeColor', markerEdgeColor);
                        end
                    end
                end
            end
        end
        
    end
    
%     plot([0 numTime/60], [27.5, 27.5], '--k')
    hold off
    xlabel('Time (hour)')
    ylabel('Contra FA-FA delay (s)')
    xlim([0 numTime/60])
    ylim([0 26])
    box off
    
    setPrint(8*2, 6*2, [plotDir, 'ContraFADelay_' fileName]);
    
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