%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.  Covariance analysis: Factor analysis -- Loading matrix using LONO 
% number of factors -- Varimax
% 
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function FACluster_v2_0(nFile)            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>    
    load([tempDatDir, fileName, '.mat'], 'dff','timePoints'); 
    load([tempDatDir, fileName, '_numFactorNoLONOActiveNeuronsSimplified.mat'], 'LONOM', 'sigNeuronsMat')
    
    numPlot           = length(LONOM);
    LMat              = cell(numPlot, 1);
    PsiMat            = cell(numPlot, 1);
    numNeuron         = size(sigNeuronsMat, 1); %#ok<NODEF>
    
    for nPlot         = 1:numPlot
        if LONOM(nPlot)  == 0
            LMat{nPlot}   = nan(numNeuron, 1);
            PsiMat{nPlot} = ones(numNeuron, 1);
        else
            LMat_nPlot    = nan(numNeuron, LONOM(nPlot));
            PsiMat_nPlot  = ones(numNeuron, LONOM(nPlot));
            
            slicedDFF     = dff(sigNeuronsMat(:, nPlot),timePoints(nPlot)+1:timePoints(nPlot)+1200); 
            slicedDFF     = bsxfun(@minus, slicedDFF, mean(slicedDFF,2));
            slicedDFF     = bsxfun(@rdivide, slicedDFF, std(slicedDFF,[],2))';
            
            [Lambda, Psi] = factoran(slicedDFF, LONOM(nPlot),'maxit',10000);
            
            LMat_nPlot(sigNeuronsMat(:, nPlot), :)   = Lambda;
            PsiMat_nPlot(sigNeuronsMat(:, nPlot)) = Psi;
            
            LMat{nPlot}   = LMat_nPlot;
            PsiMat{nPlot} = PsiMat_nPlot;
        end
    end
    
    save([tempDatDir, fileName, '_LONOLoading.mat'], 'LMat', 'PsiMat'); 
    
end