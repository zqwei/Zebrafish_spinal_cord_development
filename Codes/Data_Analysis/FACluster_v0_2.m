%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.  Covariance analysis: Factor analysis -- number of island using non
% LONO methods with selected neurons
% 
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function FACluster_v0_2(nFile)            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>    
    load([tempDatDir, fileName, '.mat'], 'dff', 'timePoints');
    load([tempDatDir, fileName, '_numFactorNoLONOActiveNeuronsSimplified.mat'], 'sigNeuronsMat'); 
    maxNumFactor      = 20;
    numPlot           = length(timePoints);
    
    numFactors        = repmat(struct('kgM', 0, 'paM', 0, 'SRMRM', 0, 'CFIM', 0), numPlot, 1);
    
    for nPlot        = 1:numPlot
        slicedDFF    = dff(:,timePoints(nPlot)+1:timePoints(nPlot)+1200); %#ok<NODEF>
        slicedDFF    = bsxfun(@minus, slicedDFF, mean(slicedDFF,2));
        slicedDFF    = bsxfun(@rdivide, slicedDFF, std(slicedDFF,[],2))';
        slicedDFF    = slicedDFF(:, sigNeuronsMat(:, nPlot)); %#ok<NODEF>     
        if ~isempty(slicedDFF)
            KG_M          = sum(eig(corr(slicedDFF))>1);
            currMaxFactor = min(KG_M, maxNumFactor);
            if currMaxFactor > 0
                numFactors(nPlot) = numFactorsWithNoCrossValidationSimplied(slicedDFF, currMaxFactor); 
            end
        end
        disp(nPlot)
    end
    
    kgM = [numFactors.kgM]; %#ok<NASGU>
    paM = [numFactors.paM]; %#ok<NASGU>
    SRMRM = [numFactors.SRMRM]; %#ok<NASGU>
    CFIM = [numFactors.CFIM]; %#ok<NASGU>
    
%     plot(timePoints(1:5:numPlot)/4/3600, kgM(1:5:numPlot), 's', ...
%         timePoints(1:5:numPlot)/4/3600, paM(1:5:numPlot), 'o', ...
%         timePoints(1:5:numPlot)/4/3600, SRMRM(1:5:numPlot), 'x', ...
%         timePoints(1:5:numPlot)/4/3600, CFIM(1:5:numPlot), '+')
%     set(gca, 'Xtick', 1:6)
%     xlim([0 timePoints(end)/4/3600])
%     xlabel('Time (hour)')
%     ylabel('# Dim')
%     legend('KG', 'PA', 'SRMR', 'CFI');
%     
%     setPrint(8, 6, [plotDir, 'numFactorActiveNeurons_', fileName], 'pdf');
    save([tempDatDir, fileName, '_numFactorNoLONOActiveNeuronsSimplified.mat'], 'kgM', 'paM', 'SRMRM', 'CFIM', '-append'); 
        
end