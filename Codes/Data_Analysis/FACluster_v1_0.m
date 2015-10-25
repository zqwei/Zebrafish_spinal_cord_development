%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.  Covariance analysis: Factor analysis -- number of island using LONO 
% methods with selected neurons
% 
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function FACluster_v1_0(nFile)            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>    
    load([tempDatDir, fileName, '.mat'], 'dff', 'timePoints');
    load([tempDatDir, fileName, '_numFactorNoLONOActiveNeuronsSimplified.mat'], 'sigNeuronsMat'); 
    
    maxNumFactor      = 12;
    numFold           = 10;  
    
    numPlot           = length(timePoints);
    EVLONO            = nan(numPlot, maxNumFactor, numFold);
    
    numSample         = 1200;
    numTest           = ceil(numSample/numFold);

    
    for nPlot        = 1:numPlot
        slicedDFF    = dff(:,timePoints(nPlot)+1:timePoints(nPlot)+1200); %#ok<NODEF>
        slicedDFF    = bsxfun(@minus, slicedDFF, mean(slicedDFF,2));
        slicedDFF    = bsxfun(@rdivide, slicedDFF, std(slicedDFF,[],2))';
        slicedDFF    = slicedDFF(:, sigNeuronsMat(:, nPlot)); %#ok<NODEF>     
        if ~isempty(slicedDFF)
            numActive     = size(slicedDFF, 2);
%             currNumFactor = min([floor(numActive/2), maxNumFactor]);
            currNumFactor = min([numActive-1, floor(numActive + 0.5 -sqrt(2*numActive + 0.25)), maxNumFactor]);
            if currNumFactor >= 1
                for nFactor         = 1:currNumFactor
                    for nFold       = 1:numFold
                        randSeqTest                  = randperm(numSample);
                        nFoldDFFTest                 = slicedDFF(randSeqTest(1:numTest), :);
                        nFoldDFFTrain                = slicedDFF(randSeqTest(numTest+1:numSample), :);
                        [lambda,psi]                 = factoran(nFoldDFFTrain, nFactor, 'scores','regression', 'rotate', 'none');
                        EVLONO(nPlot, nFactor, nFold) = LONOFA(nFoldDFFTest, lambda, psi);
                    end
                end
            end
        end
    end    
    save([tempDatDir, fileName, '_numFactorNoLONOActiveNeuronsSimplified.mat'], 'EVLONO', '-append');  
    
end