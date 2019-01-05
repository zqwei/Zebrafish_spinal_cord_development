%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.  Covariance analysis: Factor analysis -- loadin matrix evolution --
% reordering
% 
% size of factor against time
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function FACluster_v1_1_2(nFile)            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>    
    load([tempDatDir, fileName, '.mat'], 'activeNeuronMat', 'mnx'); 
    if ~exist([tempDatDir, 'EvoLoading_' fileName, '_v2.mat'], 'file'); return; end
    load([tempDatDir, 'EvoLoading_' fileName, '_v2.mat'], 'networkMat')
    numActNeuron      = sum(activeNeuronMat, 1);
    numTime           = length(numActNeuron);
    numFactor         = zeros(numTime, 1);
    factorNeuronMat   = false(size(activeNeuronMat));
    

    
    for nTime         = 1:numTime
        factorSet     = networkMat{nTime, 1};
        for nFactor   = 2:length(factorSet)
            if length(factorSet(nFactor).neuronIndex) > 1
                factorNeuronMat(factorSet(nFactor).neuronIndex, nTime) = true;
                numFactor(nTime)       = numFactor(nTime)+1;
            end
        end
    end
    
    numFactorNeuron   = sum(factorNeuronMat, 1);
    
    factorMNXTime     = sum(factorNeuronMat(~mnx, :), 1);
%     factorMNXTime     = smooth(factorMNXTime, 7);
    factorMNXTime     = find(factorMNXTime>=1, 4, 'first');
    
    factorPeakTime    = find(numFactor == max(numFactor), 4, 'first');
    
    
    if ~isempty(factorMNXTime)
        minValues = nan(length(factorMNXTime));
        minInds   = nan(length(factorMNXTime));
        for nFactor = 1:length(factorMNXTime)
            [minValues(nFactor), minInds(nFactor)] = min(abs(factorMNXTime(nFactor) - factorPeakTime));
        end
        [min_diff, ind] = min(minValues(~isnan(minValues)));
        disp([nFile, min_diff])
        factorMNXTime   = factorMNXTime(ind);
        factorPeakTime  = factorPeakTime(minInds(ind));
    end
    
    figure;
    subplot(3, 1, 1)
    hold on
    plot((1:numTime)/60, numActNeuron, 'ok')
    plot((1:numTime)/60, numFactorNeuron, 'sr')
    if ~isempty(factorMNXTime)
        gridxy(factorMNXTime/60, [], 'color', 'k', 'linestyle', '--')
    end
    gridxy(factorPeakTime/60, [], 'color', 'b', 'linestyle', '--')
    hold off
    xlabel('Time (hour)')
    ylabel('# of neuron')
    xlim([0 numTime/60])
    box off
    
    subplot(3, 1, 2)
    hold on
    plot((1:numTime)/60, numFactorNeuron./numActNeuron, 'ok')
    if ~isempty(factorMNXTime)
        gridxy(factorMNXTime/60, [], 'color', 'k', 'linestyle', '--')
    end
    gridxy(factorPeakTime/60, [], 'color', 'b', 'linestyle', '--')
    hold off
    xlabel('Time (hour)')
    ylabel('frac. neuron factored')
    xlim([0 numTime/60])
    box off
        
    subplot(3, 1, 3)
    hold on
    plot((1:numTime)/60, numFactor, 'ok')
    if ~isempty(factorMNXTime)
        gridxy(factorMNXTime/60, [], 'color', 'k', 'linestyle', '--')
    end
    gridxy(factorPeakTime/60, [], 'color', 'b', 'linestyle', '--')    
    hold off
    xlabel('Time (hour)')
    ylabel('# of Func. Com.')
    xlim([0 numTime/60])
    box off
    
    setPrint(8, 18, [plotDir, 'FASize_', fileName], 'pdf');
        
end