%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.  Covariance analysis: cross-middle line activity neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 3.3 Shuffling test build up correlation distribution for each neuron
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
%

function Data_Analysis_List3_3(nFile)
    
    addpath('../Func');
    setDir;    
    fileDirName         = fileDirNames{nFile};
    fileName            = fileNames{nFile};
    load([tempDatDir, fileName, '.mat']);

%     if ~exist([plotDir, 'Shuffle_correlation_distribution_', fileName], 'dir')
%         mkdir([plotDir, 'Shuffle_correlation_distribution_', fileName]);
%     end
        
    numShufflePoint     = 1000;
    numPlot             = 20;
    numNeuron           = size(dff,1);
    CorrDist            = zeros(numNeuron, numPlot, numShufflePoint);
    
    for nPlot                 = 1:numPlot 
        for nNeuron           = 1:numNeuron
            slicedDFF         = dff(nNeuron,timePoints(nPlot)+1:timePoints(nPlot+1));
            for nShufflePoint = 1:numShufflePoint
                testNeuron    = ceil(rand()*numNeuron);
                testPlot      = ceil(rand()*numPlot);
                testDFF       = dff(testNeuron,timePoints(testPlot)+1:timePoints(testPlot+1));
                CorrDist(nNeuron, nPlot, nShufflePoint) = corr(slicedDFF', testDFF');
            end
        end
    end
    
    
%     percentageCorr            = nan(numNeuron, numNeuron);
    mCol                      = 5;
    mRow                      = ceil(numPlot/mCol);
    pValue                    = 0.05;
    pMat                      = false(numNeuron, numNeuron);
    
    figure;
    
    for nPlot                 = 1:numPlot 
        slicedDFF             = dff(:,timePoints(nPlot)+1:timePoints(nPlot+1));   
        corrDFF               = corr(slicedDFF'); 
        for nNeuron           = 1:numNeuron
            distNeuron        = squeeze(CorrDist(nNeuron, nPlot,:));
            for mNeuron       = 1:numNeuron
%                 percentageCorr(nNeuron, mNeuron) = sum(corrDFF(nNeuron, mNeuron)>distNeuron)/numShufflePoint;
                percentageCorr         = sum(corrDFF(nNeuron, mNeuron)>distNeuron)/numShufflePoint;
                pMat(nNeuron, mNeuron) = min(percentageCorr, 1-percentageCorr) < pValue/2;
            end
        end
        subplot(mRow, mCol, nPlot)
%         imagesc(percentageCorr, [0 1]);
        imagesc(pMat);
        colormap(gray)
    end    
    
    setPrint(8*mCol, 6*mRow, [plotDir, 'Shuffle_correlation_distribution_', fileName], 'pdf');
    save([tempDatDir, fileName, '_Shuffle_correlation_distribution_PMat.mat'], 'pMat');

end