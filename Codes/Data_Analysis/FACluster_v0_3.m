%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.  Covariance analysis: Factor analysis -- number of island using non
% LONO methods with active neurons
%
%
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
%


function FACluster_v0_3(nFile)
    addpath('../Func');
    setDir;
    fileName          = fileNames{nFile}; %#ok<USENS>
    load([tempDatDir, fileName, '.mat'], 'dff','timePoints','activeNeuronMat', 'timeStep'); 
    load([tempDatDir, 'FALONO_', fileName, '.mat'], 'uncorrectedLONOM');

    LONOM             = uncorrectedLONOM;
    numPlot           = length(LONOM);
    LMat              = cell(numPlot, 1);
    PsiMat            = cell(numPlot, 1);
    numNeuron         = size(activeNeuronMat, 1); %#ok<NODEF>

    for nPlot         = 1:numPlot
        if LONOM(nPlot)  == 0
            LMat{nPlot}   = nan(numNeuron, 1);
            PsiMat{nPlot} = ones(numNeuron, 1);
        else
            LMat_nPlot    = nan(numNeuron, LONOM(nPlot));
            PsiMat_nPlot  = ones(numNeuron, 1);

            slicedDFF     = dff(activeNeuronMat(:, nPlot),timePoints(nPlot)+1:timePoints(nPlot)+timeStep);
            % remove data with twitch times
            slicedDFF(:, sum(isnan(slicedDFF))>0)     = [];
            slicedDFF     = bsxfun(@minus, slicedDFF, mean(slicedDFF,2));
            slicedDFF     = bsxfun(@rdivide, slicedDFF, std(slicedDFF,[],2))';

            [Lambda, Psi] = factoran(slicedDFF, LONOM(nPlot),'maxit',10000);

            LMat_nPlot(activeNeuronMat(:, nPlot), :)   = Lambda;
            PsiMat_nPlot(activeNeuronMat(:, nPlot)) = Psi;

            LMat{nPlot}   = LMat_nPlot;
            PsiMat{nPlot} = PsiMat_nPlot;
        end
    end

    save([tempDatDir, 'LONOLoading_' fileName, '.mat'], 'LMat', 'PsiMat');
end
