%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5.  MNX neurons -- active neurons
% 
%  EV and activation as a function of time
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function MNX_v0_6(nFile)            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>    
    load([tempDatDir, fileName, '.mat'], 'mnx', 'tracks', 'side', 'leafOrder', 'slicedIndex', 'activeNeuronMat'); 
    load([tempDatDir, 'EvoLoading_' fileName, '_v2.mat'], 'networkMat')
    
    if ~exist('mnx', 'var')
        return;
    end    
        
    numTime       = size(networkMat, 1); %#ok<*NODEF>
    meanDelayTime = nan(numTime, 1);
    stdDelayTime  = nan(numTime, 1);
    
    
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
        delaySet        = [];
        for nFactor     = find(FactorIndex)
            for mFactor = find(FactorIndex)
                if nFactor < mFactor
                    delayTime = networkCorr.delayMat(nFactor-1, mFactor-1);
                    delayTime = abs(delayTime);
                    isContra  = networkCorr.IpsiIndex(nFactor-1, mFactor-1) == 0;
                    if isContra
                        delaySet  = [delayTime; delaySet];
                    end
                end
            end
        end
        
        delaySet(isnan(delaySet))  = 50;
        meanDelayTime(nTime)       = mean(delaySet);
        stdDelayTime(nTime)        = std(delaySet);
    end
    
    figure; hold on;
    
    fracMNXPos    = mean(activeNeuronMat(mnx==1, :));
    fracMNXNeg    = mean(activeNeuronMat(mnx==0, :));
    
    yyaxis right
    plot((1:numTime)/60, smooth(fracMNXPos/max(fracMNXPos), 11, 'rlowess'), '-b', 'linewid', 2)
    if max(fracMNXNeg)>0
        plot((1:numTime)/60, smooth(fracMNXNeg/max(fracMNXNeg), 11, 'rlowess'), '-r', 'linewid', 2)
    end
    ylim([0 1])
    ylabel('Normalized fraction neuron')
    
    yyaxis left
%     plot((1:numTime)/60, smooth(meanDelayTime, 11, 'rlowess'), '-k', 'linewid', 1)     
%     plot((1:numTime)/60, smooth(meanDelayTime-stdDelayTime, 11, 'rlowess'), '-k', 'linewid', 0.5)
%     plot((1:numTime)/60, smooth(meanDelayTime+stdDelayTime, 11, 'rlowess'), '-k', 'linewid', 0.5)
%     plot((1:numTime)/60, stdDelayTime, '-k', 'linewid', 1)
%     shadedErrorBar((1:numTime)/60, smooth(meanDelayTime, 11), smooth(stdDelayTime, 11), {'-k', 'linewid', 0.5}, 0.5)
    shadedErrorBar((1:numTime)/60, meanDelayTime, stdDelayTime, {'-k', 'linewid', 1}, 0.5)
    
    hold off
    xlabel('Time (hour)')
    ylabel('Contra FA-FA delay (s)')
    xlim([0 numTime/60])
    ylim([0 35])
    box off
    
    setPrint(8, 6, [plotDir, 'ContraFADelayFracMNXNeuron_' fileName]);
    
end