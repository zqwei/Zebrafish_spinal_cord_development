%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.  Covariance analysis: Factor analysis -- number of island using non
% LONO methods with all neurons
% 
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function FACluster_v0_0(nFile)            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>    
    load([tempDatDir, fileName, '.mat'], 'dff', 'timePoints');
    maxNumFactor      = 20;
    numPlot           = length(timePoints);
    
    numFactors        = repmat(struct('kgM', 0, 'paM', 0, 'SRMRM', 0, 'CFIM', 0), numPlot, 1);
    
    for nPlot        = 1:5:numPlot
        slicedDFF    = dff(:,timePoints(nPlot)+1:timePoints(nPlot)+1200); %#ok<NODEF>
        slicedDFF    = bsxfun(@minus, slicedDFF, mean(slicedDFF,2));
        slicedDFF    = bsxfun(@rdivide, slicedDFF, std(slicedDFF,[],2))';
        numFactors(nPlot) = numFactorsWithNoCrossValidationSimplied(slicedDFF, maxNumFactor); 
        disp(nPlot)
    end
    
    kgM = [numFactors.kgM];
    paM = [numFactors.paM];
    SRMRM = [numFactors.SRMRM];
    CFIM = [numFactors.CFIM];
    
    plot(timePoints(1:5:numPlot)/4/3600, kgM(1:5:numPlot), 'o', ...
        timePoints(1:5:numPlot)/4/3600, paM(1:5:numPlot), 'o', ...
        timePoints(1:5:numPlot)/4/3600, SRMRM(1:5:numPlot), 'o', ...
        timePoints(1:5:numPlot)/4/3600, CFIM(1:5:numPlot), 'o')
    xlim([0 timePoints(end)/4/3600])
    xlabel('Time (hour)')
    ylabel('# Dim')
    legend('KG', 'PA', 'SRMR', 'CFI');
    
    setPrint(8, 6, [plotDir, 'numFactorAllNeurons_', fileName], 'pdf');
        
end