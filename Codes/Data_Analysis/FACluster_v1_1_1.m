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


function FACluster_v1_1_1(nFile)            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>    
    load([tempDatDir, 'EvoLoading_' fileName, '_v2.mat'], 'networkMat')
    numTime           = length(networkMat);
    numFactorSize     = nan(numTime, 2);
    
    
    for nTime         = 1:numTime
        factorSet     = networkMat{nTime, 1};
        leftSize      = [];
        rightSize     = [];
        for nFactor   = 2:length(factorSet)
            if length(factorSet(nFactor).neuronIndex) > 1
                if factorSet(nFactor).y>0
                    leftSize = [leftSize; length(factorSet(nFactor).neuronIndex)];
                else
                    rightSize = [rightSize; length(factorSet(nFactor).neuronIndex)];
                end
            end
        end
        if ~isempty(leftSize)
            numFactorSize(nTime, 1) = max(leftSize);
        end
        if ~isempty(rightSize)
            numFactorSize(nTime, 2) = max(rightSize);
        end
    end
    
    figure;
    plot((1:numTime)/60, numFactorSize, 'o')
    xlabel('Time (hour)')
    ylabel('Size of Func. Com.')
    xlim([0 numTime/60])
    box off
    
    setPrint(8, 6, [plotDir, 'maxFASizeTime_', fileName], 'pdf');
end