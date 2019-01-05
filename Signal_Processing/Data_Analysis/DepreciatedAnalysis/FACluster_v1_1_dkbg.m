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


function FACluster_v1_1_dkbg(nFile)            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>    
    load([tempDatDir, fileName, '.mat'], 'activeNeuronMat'); 
    load([tempDatDir, 'EvoLoading_' fileName, '.mat'], 'networkMat')
    numActNeuron      = sum(activeNeuronMat, 1);
    numTime           = length(numActNeuron);
    numFactor         = zeros(numTime, 1);
    numFactorSize     = nan(numTime, 10);
    numFactorNeuron   = numActNeuron;
    
    
    for nTime         = 1:numTime
        factorSet     = networkMat{nTime};
        for nFactor   = 2:length(factorSet)-1
            if length(factorSet{nFactor}.neuronIndex) == 1
                numFactorNeuron(nTime) = numFactorNeuron(nTime) - 1;
            else
                numFactor(nTime)       = numFactor(nTime)+1;
                numFactorSize(nTime, numFactor(nTime)) = length(factorSet{nFactor}.neuronIndex);
            end
        end
    end
    
    figure;
    plot((1:numTime)/60, numFactorSize, 'ow')
    xlabel('Time (hour)')
    ylabel('Size of Func. Com.')
    xlim([0 numTime/60])
    box off
    
%     subplot(3, 1, 2)
%     hold on
%     plot((1:numTime)/60, numActNeuron, 'ow')
%     plot((1:numTime)/60, numFactorNeuron, 'sr')
%     hold off
%     xlabel('Time (hour)')
%     ylabel('# of neuron')
%     xlim([0 numTime/60])
%     box off
%     
%     subplot(3, 1, 3)
%     plot((1:numTime)/60, numFactor, 'ow')
%     xlabel('Time (hour)')
%     ylabel('# of Func. Com.')
%     xlim([0 numTime/60])
%     box off
%     
    setPrint(8, 6, [plotDir, 'FASize_', fileName '_dkbg'], 'pdf');
end