%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  Covariance analysis: factor analysis
%     -- using PSDPeakTime to filered out the non-active units and low
%     frequency units
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

function Data_Analysis_List2_0_9(nFile)
    % load data
    addpath('../Func');
    setDir;
    
    fileDirName       = fileDirNames{nFile}; %#ok<USENS,NASGU>
    fileName          = fileNames{nFile}; %#ok<USENS>
    
    load([tempDatDir, fileName, '.mat']); %dff, timePoints
    load([tempDatDir, fileName, '_PSDPeakTime.mat'], 'PSDPeakTime')
    load([tempDatDir, fileName, '_numFactorNoLONOActiveNeuronsSimplified.mat'], 'numFactors', 'nActiveUnit');
    
    maxNumFactor  = 10;

    numFold       = 10;
    numPlot       = length(timePoints)-1;
    matEV         = nan(numPlot, maxNumFactor, numFold);
    numUnit       = size(dff,1);
    matEVSingleUnit = nan(numPlot, maxNumFactor, numUnit, numFold);



    for nPlot        = numPlot:-1:1
        activeIndex  = PSDPeakTime(:, nPlot)>0 & ~isnan(PSDPeakTime(:, nPlot)); %#ok<NODEF>
        numActive    = sum(activeIndex);
        slicedDFF    = dff(activeIndex, timePoints(nPlot)+1:timePoints(nPlot+1)); % only take the active units
        slicedDFF    = bsxfun(@minus, slicedDFF, mean(slicedDFF,2));
        slicedDFF    = bsxfun(@rdivide, slicedDFF, std(slicedDFF,[],2));
        numSample    = size(slicedDFF, 2);
        numTest      = ceil(numSample * 0.1);
        currFold     = min([floor(numActive/2), numFold]);
        currNumFactor = min([numActive, floor(numActive + 0.5 -sqrt(2*numActive + 0.25)), maxNumFactor]);
        for nFactor         = 1:currNumFactor
            for nFold       = 1:currFold
                randSeqTest                  = randperm(numSample);
                nFoldDFFTest                 = slicedDFF(:,randSeqTest(1:numTest))';
                nFoldDFFTrain                = slicedDFF(:,randSeqTest(numTest+1:numSample))';
                [lambda,psi]                 = factoran(nFoldDFFTrain, nFactor, 'scores','regression', 'rotate', 'none');
                matEV(nPlot, nFactor, nFold) = LONOFA(nFoldDFFTest, lambda, psi);
                matEVSingleUnit(nPlot, nFactor, activeIndex, nFold) = LONOFASingleUnitEV(nFoldDFFTest, lambda, psi);
                display([nPlot, nFactor, nFold])
            end
        end
    end

    % figure;
    % plot(1:maxNumFactor, mean(matEV,3))
    % errorbar(repmat(timePoints(1:end-1),maxNumFactor,1), mean(matEV,3)', std(matEV,[],3)');

    save([tempDatDir, fileName, '_PSDPeakTimeFAEV.mat'], 'matEV');
% 
%     max_ev_units = squeeze(max(squeeze(mean(matEVSingleUnit,4)),[],2));
%     figure ;
%     for ii = 1:3
%         subplot(1,3,ii);
%         set(0,'DefaultAxesColorOrder',distinguishable_colors(20))
%         plot(timePoints(1:end-1)/4/3600,max_ev_units(:,idx==ii)),ylim([0 1]) %#ok<COLND>
%     end
    save([tempDatDir, fileName, '_PSDPeakTimeFAEV.mat'], 'matEVSingleUnit','-append');
    
end