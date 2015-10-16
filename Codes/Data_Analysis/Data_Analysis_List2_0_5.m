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

function Data_Analysis_List2_0_5(nFile, dff_threshould)
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>
    
    load([tempDatDir, fileName, '.mat']);
    load([tempDatDir, fileName, '_PSDPeakTime.mat'], 'PSDPeakTime')
    maxNumFactor      = 20;
    
    % Computing numFactors from non-LONO methods
    numPlot           = length(timePoints)-1;
    numFactors        = repmat(struct('kgM', 0, 'varM', 0, 'vMapM', 0, 'paM', 0, ...
                        'pChisqTestM', 0, 'AICM', 0, 'BICM', 0, 'CAICM', 0, ...
                        'SRMRM', 0, 'RMSEAM', 0, 'CFIM', 0), numPlot, 1);
                    
%     dff_threshould   = 0.03;
    nSkipTime        = true(numPlot, 1);
    nActiveUnit      = zeros(numPlot, 1);
    
    for nPlot        = 1:numPlot
        slicedDFF    = dff(:,timePoints(nPlot)+1:timePoints(nPlot+1)); %#ok<NODEF>
        activeIndex  = std(slicedDFF, [], 2) > dff_threshould;
        slicedDFF    = slicedDFF(activeIndex, :);
        slicedDFF    = slicedDFF(activeIndex, :);
        slicedDFF    = bsxfun(@minus, slicedDFF, mean(slicedDFF,2));
        if sum(activeIndex) > 3
            nSkipTime(nPlot)  = false;
            nActiveUnit(nPlot) = sum(activeIndex);
            numFactors(nPlot) = numFactorsWithNoCrossValidation(slicedDFF', maxNumFactor);       
        end
    end
    
    save([tempDatDir, fileName, '_numFactorNoLONOActiveNeuronsThres.mat'], 'numFactors', 'nSkipTime', 'nActiveUnit'); 
    
end