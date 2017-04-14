%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5.  MNX neurons -- active neurons
% 
%  EV and activation as a function of time
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function MNX_v0_4(nFile)            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>    
    mnx               = [];
    load([tempDatDir, fileName, '.mat'], 'mnx', 'tracks', 'side', 'leafOrder', 'slicedIndex', 'activeNeuronMat'); 
    load([tempDatDir, 'EV_' fileName, '.mat'], 'EVMat')
    load([tempDatDir, 'EvoLoading_' fileName, '_v2.mat'], 'networkMat')
    
    if isempty(mnx) || sum(~mnx) == 0
        return;
    end
        
    numTime           = size(EVMat, 2); %#ok<*NODEF>
    numNeuron         = size(EVMat, 1);     
    mnxNegIndex       = find(~mnx);
    numMNXNeuron      = length(mnxNegIndex);
    
    neuronFactor      = nan(numMNXNeuron, numTime, 2); % first col: corr, second col: delay
    factorIpsiFactor  = nan(numMNXNeuron, numTime, 10, 2); % maximum 10 factor ipsi
    factorContraFactor= nan(numMNXNeuron, numTime, 10, 2); % maximum 10 factor contra
    
    timeLag           = 5;
    
    for mNeuron       = 1:numMNXNeuron
        nNeuron       = mnxNegIndex(mNeuron);
        
        % current time
        
        for nTime     = 1:numTime
            networkTime    = networkMat{nTime, 1}; %#ok<*USENS>
            networkCorr    = networkMat{nTime, 2};
            numUnitsFactor = [networkTime.nUnit];
            FactorIndex    = numUnitsFactor > 1;
            FactorIndex(1) = false;
            
            factorNeuron   = [];
            for nFactor    = find(FactorIndex)
                if sum(networkTime(nFactor).neuronIndex == nNeuron) > 0
                    factorNeuron = nFactor;
                    % disp([nNeuron, nTime, nFactor])
                    currFactorNeuron   = networkTime(nFactor).neuronIndex;
                end
            end
            
            % if the neuron is not belong to any factor, break the
            % computation
            if isempty(factorNeuron)
                continue;
            end
            
            neuronFactor(mNeuron, nTime, :)   = networkTime(factorNeuron).neuronCCMat(nNeuron, [2 1]);
            
            FactorIndex(factorNeuron) = false;
            nPairContra      = 0;
            nPairIpsi        = 0;
            for mFactor      = find(FactorIndex)    
                if mFactor > factorNeuron
                    if networkCorr.IpsiIndex(factorNeuron-1, mFactor-1) == 1
                        nPairIpsi           = nPairIpsi + 1;
                        factorIpsiFactor(mNeuron, nTime, nPairIpsi, 1) = networkCorr.corrMat(factorNeuron-1, mFactor-1);
                        factorIpsiFactor(mNeuron, nTime, nPairIpsi, 2) = networkCorr.delayMat(factorNeuron-1, mFactor-1);
                    elseif networkCorr.IpsiIndex(factorNeuron-1, mFactor-1) == 0
                        nPairContra         = nPairContra + 1;
                        factorContraFactor(mNeuron, nTime, nPairContra, 1) = networkCorr.corrMat(factorNeuron-1, mFactor-1);
                        factorContraFactor(mNeuron, nTime, nPairContra, 2) = networkCorr.delayMat(factorNeuron-1, mFactor-1);
                    end
                else
                    if networkCorr.IpsiIndex(mFactor-1, factorNeuron-1) == 1
                        nPairIpsi           = nPairIpsi + 1;
                        factorIpsiFactor(mNeuron, nTime, nPairIpsi, 1) = networkCorr.corrMat(mFactor-1, factorNeuron-1);
                        factorIpsiFactor(mNeuron, nTime, nPairIpsi, 2) = networkCorr.delayMat(mFactor-1, factorNeuron-1);
                    elseif networkCorr.IpsiIndex(mFactor-1, factorNeuron-1) == 0
                        nPairContra         = nPairContra + 1;
                        factorContraFactor(mNeuron, nTime, nPairContra, 1) = networkCorr.corrMat(mFactor-1, factorNeuron-1);
                        factorContraFactor(mNeuron, nTime, nPairContra, 2) = networkCorr.delayMat(mFactor-1, factorNeuron-1);
                    end
                end
            end
            
            % previous time point in time lag
            for mTime    = nTime - 1:-1:nTime - timeLag
                if mTime > 0 && isnan(neuronFactor(mNeuron, mTime, 1)) % if mTime is never computed
                    networkTime    = networkMat{mTime, 1}; %#ok<*USENS>
                    networkCorr    = networkMat{mTime, 2};
                    numUnitsFactor = [networkTime.nUnit];
                    FactorIndex    = numUnitsFactor > 1;
                    FactorIndex(1) = false;
                    
                    numOverlap     = nan(length(numUnitsFactor), 1);
                    for nFactor    = find(FactorIndex)
                        numOverlap(nFactor) = length(intersect(networkTime(nFactor).neuronIndex, currFactorNeuron));
                    end
                    
                    [maxOverlap, maxOverlapInd] = max(numOverlap);
                    
                    % if the neuron is not belong to any factor, break the
                    % computation
                    if isnan(maxOverlap)
                        continue;
                    end
                    
                    factorNeuron = maxOverlapInd(1);

                    neuronFactor(mNeuron, mTime, :)   = networkTime(factorNeuron).neuronCCMat(nNeuron, [2 1]);

                    FactorIndex(factorNeuron) = false;
                    nPairContra      = 0;
                    nPairIpsi        = 0;
                    for mFactor      = find(FactorIndex)    
                        if mFactor > factorNeuron
                            if networkCorr.IpsiIndex(factorNeuron-1, mFactor-1) == 1
                                nPairIpsi           = nPairIpsi + 1;
                                factorIpsiFactor(mNeuron, mTime, nPairIpsi, 1) = networkCorr.corrMat(factorNeuron-1, mFactor-1);
                                factorIpsiFactor(mNeuron, mTime, nPairIpsi, 2) = networkCorr.delayMat(factorNeuron-1, mFactor-1);
                            elseif networkCorr.IpsiIndex(factorNeuron-1, mFactor-1) == 0
                                nPairContra         = nPairContra + 1;
                                factorContraFactor(mNeuron, mTime, nPairContra, 1) = networkCorr.corrMat(factorNeuron-1, mFactor-1);
                                factorContraFactor(mNeuron, mTime, nPairContra, 2) = networkCorr.delayMat(factorNeuron-1, mFactor-1);
                            end
                        else
                            if networkCorr.IpsiIndex(mFactor-1, factorNeuron-1) == 1
                                nPairIpsi           = nPairIpsi + 1;
                                factorIpsiFactor(mNeuron, mTime, nPairIpsi, 1) = networkCorr.corrMat(mFactor-1, factorNeuron-1);
                                factorIpsiFactor(mNeuron, mTime, nPairIpsi, 2) = networkCorr.delayMat(mFactor-1, factorNeuron-1);
                            elseif networkCorr.IpsiIndex(mFactor-1, factorNeuron-1) == 0
                                nPairContra         = nPairContra + 1;
                                factorContraFactor(mNeuron, mTime, nPairContra, 1) = networkCorr.corrMat(mFactor-1, factorNeuron-1);
                                factorContraFactor(mNeuron, mTime, nPairContra, 2) = networkCorr.delayMat(mFactor-1, factorNeuron-1);
                            end
                        end
                    end
                else
                    break;
                end
            end
        end
        
    end
    
    
    plotFunc(squeeze(neuronFactor(:, :, 1)), [plotDir, 'SigFitEVMNXNegActiveNeuron_', fileName '_neuron_FA_corr'], 1);
    plotFunc(squeeze(neuronFactor(:, :, 2)), [plotDir, 'SigFitEVMNXNegActiveNeuron_', fileName '_neuron_FA_delay'], 1);
%     plotFunc(squeeze(factorIpsiFactor(:, :, :, 1)), [plotDir, 'SigFitEVMNXNegActiveNeuron_', fileName '_FA_IpsiFA_corr'], yLimMax);
%     plotFunc(squeeze(factorIpsiFactor(:, :, :, 2)), [plotDir, 'SigFitEVMNXNegActiveNeuron_', fileName '_FA_IpsiFA_delay'], yLimMax);
    plotFunc(squeeze(factorContraFactor(:, :, :, 1)), [plotDir, 'SigFitEVMNXNegActiveNeuron_', fileName '_FA_ContraFA_corr'], 1);
    plotFunc(squeeze(factorContraFactor(:, :, :, 2)), [plotDir, 'SigFitEVMNXNegActiveNeuron_', fileName '_FA_ContraFA_delay'], 5);

    
    
end


function plotFunc(plotMat, fileName, yLimMax)
    
    numNeuron   = size(plotMat, 1);
    m           = ceil(sqrt(numNeuron));
    figure;
    numTime     = size(plotMat, 2);
    timeIndex   = 1:numTime;
    
    for nNeuron       = 1:numNeuron
        subplot(m, m, nNeuron)
        if ismatrix(plotMat)
            plotMatNeuron  = plotMat(nNeuron, :); 
        else
            plotMatNeuron  = squeeze(plotMat(nNeuron, :, :)); 
        end
        plot(timeIndex/60, abs(plotMatNeuron), '.k');
        xlim([0 numTime/60]) 
        ylim([0 yLimMax])
        xlabel('Time (hour)')
        box off
    end
    setPrint(8*m, 6*m, fileName, 'pdf');
end