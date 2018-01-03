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


function FACluster_v0_4(nFile)
    addpath('../Func');
    setDir;
    fileName          = fileNames{nFile}; %#ok<USENS>
    load([tempDatDir, fileName, '.mat'], 'sideSplitter', 'side', 'tracks', 'timePoints', 'activeNeuronMat', 'dff', 'timePoints', 'timeStep'); 
    load([tempDatDir, 'LONOLoading_' fileName, '.mat'], 'LMat', 'PsiMat')
    numTime           = length(LMat); %#ok<NODEF>
    CorrectedLMat     = LMat;
    pValue            = 0.05;
    EVThres           = 0.05;

    for nTime         = 1:numTime
        LMat              = CorrectedLMat{nTime}; % threshold by 0.3
        nanIndex          = ~isnan(LMat(:,1));
        nanLMat           = LMat(nanIndex, :);
        psi               = PsiMat{nTime};
        psi               = psi(nanIndex, :);

        if sum(nanIndex)>0
            slicedDFF         = dff(:,timePoints(nTime)+1:timePoints(nTime)+timeStep); %#ok<NODEF>
            % remove data with twitch times
            slicedDFF(:, sum(isnan(slicedDFF))>0)     = [];
            slicedDFF         = bsxfun(@minus, slicedDFF, mean(slicedDFF,2));
            slicedDFF         = bsxfun(@rdivide, slicedDFF, std(slicedDFF,[],2))';
            slicedDFF         = slicedDFF(:, nanIndex);

            [~, corrSig]      = corr(slicedDFF);
            nanLMat(sum(corrSig<pValue)==0, :) = 0;
            psi(sum(corrSig<pValue)==0)  = 1;

            EVSingleUnit      = LONOFASingleUnitSingleFactorEV(slicedDFF, nanLMat, psi);
            nanLMat(EVSingleUnit<EVThres) = 0;
            LMat(nanIndex, :)    = nanLMat;
            CorrectedLMat{nTime} = LMat;
        end

    end

    save([tempDatDir, 'LONOLoading_' fileName, '.mat'], 'CorrectedLMat', '-append');
end
