%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.  Covariance analysis: Factor analysis -- Loading matrix using LONO 
% number of factors -- Varimax -- correction of loading matrix -- comparing
% number of sig. Neurons vs neurons in factors (they should be the same)
% 
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function FACluster_v2_3(nFile)            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>    
    load([tempDatDir, fileName, '.mat'], 'timePoints', 'dff'); 
    load([tempDatDir, fileName, '_numFactorNoLONOActiveNeuronsSimplified.mat'], 'sigNeuronsMat', 'LONOM')
    load([tempDatDir, fileName, '_LONOLoading.mat'], 'CorrectedLMat', 'numDiffNeurons');
    
    numPlot           = length(LONOM);
    
    for nPlot         = 1:numPlot    
        if numDiffNeurons(nPlot)>0
            LMatnTime     = CorrectedLMat{nPlot};
            FANeurons     = sum(LMatnTime>0, 2)>0;
            missingNeuron = (sigNeuronsMat(:, nPlot)) & (~FANeurons);
                        
            if sum(missingNeuron) == 2
                slicedDFF     = dff(missingNeuron,timePoints(nPlot)+1:timePoints(nPlot)+1200); 
                slicedDFF     = bsxfun(@minus, slicedDFF, mean(slicedDFF,2));
                slicedDFF     = bsxfun(@rdivide, slicedDFF, std(slicedDFF,[],2))';
                slicedDFF     = [slicedDFF, randn(length(slicedDFF), 1)];
                [Lambda, Psi] = factoran(slicedDFF, 1,'maxit',10000);
                LMatnTime(missingNeuron, LONOM(nPlot)+1) = Lambda(1:2);
                LONOM(nPlot)       = LONOM(nPlot)+1;
            end
            CorrectedLMat{nPlot} = LMatnTime;
        end
    end
    
    save([tempDatDir, fileName, '_LONOLoading.mat'], 'CorrectedLMat', '-append')
    save([tempDatDir, fileName, '_numFactorNoLONOActiveNeuronsSimplified.mat'], 'LONOM', '-append')

end