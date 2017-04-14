%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5.  MNX neurons -- active neurons
% 
%  EV and activation as a function of time
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function MNX_v0_5(nFile)            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>    
    mnx               = [];
    load([tempDatDir, fileName, '.mat'], 'mnx', 'tracks', 'side', 'leafOrder', 'slicedIndex', 'activeNeuronMat'); 
    load([tempDatDir, 'EvoLoading_' fileName, '_v2.mat'], 'networkMat')
    
    if isempty(mnx) || sum(~mnx) == 0
        return;
    end
    
    colorSet      = [0.8 0.8 0.8; 0.2 0.2 1.0; 0.0 0.0 0.0];
        
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
                    if isContra
                        markerColor = colorSet((FactorMNX(nFactor) < 1) + (FactorMNX(mFactor) < 1) + 1, :);
                        if isnan(delayTime)
                            plot(nTime/60, 28+rand(), '.', 'color', markerColor);
                        else
                            plot(nTime/60, delayTime, '.', 'color', markerColor);
                        end
                    end
                end
            end
        end
        
    end
    
    plot([0 numTime/60], [27.5, 27.5], '--k')
    hold off
    xlabel('Time (hour)')
    ylabel('Contra FA-FA delay (s)')
    xlim([0 numTime/60])
    ylim([0 30])
    box off
    
    setPrint(8, 6, [plotDir, 'ContraFADelay_' fileName]);
    
end