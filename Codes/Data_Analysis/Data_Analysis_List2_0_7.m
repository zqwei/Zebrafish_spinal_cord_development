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

function Data_Analysis_List2_0_7(nFile)
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>
    
    load([tempDatDir, fileName, '.mat']);
    load([tempDatDir, fileName, '_PSDPeakTime.mat'], 'PSDPeakTime')
    maxNumFactor      = 15;
    
    % Computing numFactors from non-LONO methods
    numPlot           = length(timePoints)-1;
                    
    numFactors        = repmat(struct('kgM', 0, 'SRMRM', 0, 'CFIM', 0), numPlot, 1);
                    
    nSkipTime        = true(numPlot, 1);
    nActiveUnit      = zeros(numPlot, 1);
    
    for nPlot        = 1:numPlot
        slicedDFF    = dff(:,timePoints(nPlot)+1:timePoints(nPlot+1)); %#ok<NODEF>
        activeIndex  = PSDPeakTime(:, nPlot)>0 & ~isnan(PSDPeakTime(:, nPlot)); %#ok<NODEF>
        slicedDFF    = slicedDFF(activeIndex, :);
        slicedDFF    = bsxfun(@minus, slicedDFF, mean(slicedDFF,2));
        slicedDFF    = bsxfun(@rdivide, slicedDFF, std(slicedDFF,[],2))';
        if sum(activeIndex) > 3
            nSkipTime(nPlot)  = false;
            nActiveUnit(nPlot) = sum(activeIndex);
            numFactors(nPlot) = numFactorsWithNoCrossValidationSimplied(slicedDFF, maxNumFactor);            
        end
    end
    
    save([tempDatDir, fileName, '_numFactorNoLONOActiveNeuronsSimplified.mat'], 'numFactors', 'nSkipTime', 'nActiveUnit'); 
    
end