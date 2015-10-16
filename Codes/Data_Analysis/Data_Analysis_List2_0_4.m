%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  Covariance analysis: factor analysis
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.0.4 number of factors using non LONO methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

function Data_Analysis_List2_0_4(nFile)
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>
    
    load([tempDatDir, fileName, '.mat']);
    maxNumFactor      = 30; % 45 for first dataset
    
    % Computing numFactors from non-LONO methods
    numPlot           = length(timePoints)-1;
    numFactors        = repmat(struct('kgM', 0, 'varM', 0, 'vMapM', 0, 'paM', 0, ...
                        'pChisqTestM', 0, 'AICM', 0, 'BICM', 0, 'CAICM', 0, ...
                        'SRMRM', 0, 'RMSEAM', 0, 'CFIM', 0), numPlot, 1);
    
    for nPlot        = 1:numPlot
        disp(nPlot)
        slicedDFF    = dff(:,timePoints(nPlot)+1:timePoints(nPlot+1)); %#ok<NODEF>
        slicedDFF    = slicedDFF(activeIndex, :);
        slicedDFF    = bsxfun(@minus, slicedDFF, mean(slicedDFF,2));
        slicedDFF    = bsxfun(@rdivide, slicedDFF, std(slicedDFF,[],2));
        numFactors(nPlot) = numFactorsWithNoCrossValidation(slicedDFF', maxNumFactor);
    end
    
    save([tempDatDir, fileName, '_numFactorNoLONO.mat'], 'numFactors'); 
    
end