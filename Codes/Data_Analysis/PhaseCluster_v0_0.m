%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7.  Phase
% 
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
%

function PhaseCluster_v0_0(nFile, nTime)

    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>    
    load([tempDatDir, fileName, '.mat'], 'dff', 'side', 'timePoints'); 
    load([tempDatDir, fileName, '_LONOLoading.mat'], 'CorrectedLMat', 'newPsiMat')   
    load([tempDatDir, fileName, '_numFactorNoLONOActiveNeuronsSimplified.mat'], 'activeNeuronMat'); 
    numTime           = length(CorrectedLMat); %#ok<USENS>
    corrThres         = 0.3;
    
%     nTime         = 120;
    slicedDFF     = dff(:, timePoints(nTime)+1:timePoints(nTime)+1200); 
    slicedDFF     = bsxfun(@minus, slicedDFF, mean(slicedDFF,2));
    slicedDFF     = bsxfun(@rdivide, slicedDFF, std(slicedDFF,[],2))';
    LMat          = CorrectedLMat{nTime};
    Ph            = newPsiMat{nTime}(:, 1);
    LMat(LMat<corrThres) = 0;
    LMat(isnan(LMat)) = 0;
    % remove non-dominate factors
    [~, maxFactorPerNeuronIndex] = max(LMat(sum(LMat, 2)>0, :), [], 2);
    sideRemoveList  = histc(maxFactorPerNeuronIndex, 1:size(LMat, 2)) <2;    
    LMat(:, sideRemoveList) = [];
    sideRemoveList = [];
    for nFactor = 1:size(LMat, 2)
        neuronFactor = LMat(:, nFactor)>0; 
        if length(unique(side(neuronFactor)))>1
            sideRemoveList   = [sideRemoveList; nFactor]; %#ok<AGROW>
        end
    end
    LMat(:, sideRemoveList) = [];
    estX        = estFactor(LMat, Ph, slicedDFF);
    
    maxLags     = 40;
    numStd      = 3;
    fs          = 4;
    
    peakLagMat  = nan(size(slicedDFF, 2), size(estX, 1));
    peakCorrMat = nan(size(slicedDFF, 2), size(estX, 1));
    
    for nFactor = 1:size(estX, 1)
        for nNeuron = 1:size(slicedDFF, 2)
            [xcf, lags, bounds] = crosscorr(estX(nFactor,:), slicedDFF(:,nNeuron), maxLags, numStd);
            [pks, locs]         = findpeaks(xcf, lags, 'MinPeakHeight', bounds(1), 'MinPeakProminence', 0.1,'SortStr','descend');
            if ~isempty(locs)
                peakLagMat(nNeuron, nFactor)  = locs(1)/fs;
                peakCorrMat(nNeuron, nFactor) = pks(1);
            end
        end
    end
    
    figure; 
    imagesc(peakLagMat)
    xlim([0.5 size(estX, 1)+0.5])
    ylim([0.5 size(slicedDFF, 2)+0.5])
    set(gca, 'Ytick', 1:size(slicedDFF, 2))
    
    
    peakClusterLag  = nan(size(estX, 1), size(estX, 1));
    for nFactor = 1:size(estX, 1)
        for mFactor = 1:size(estX, 1)
            [xcf, lags, bounds] = crosscorr(estX(mFactor,:), estX(nFactor,:), maxLags, numStd);
            [pks, locs]         = findpeaks(xcf, lags, 'MinPeakHeight', bounds(1), 'MinPeakProminence', 0.1,'SortStr','descend');
            if ~isempty(locs)
                peakClusterLag(mFactor, nFactor)  = locs(1)/fs;
%                 peakCorrMat(mFactor, nFactor) = pks(1);
            end
        end
    end
    
    figure;
    imagesc(LMat)
    
    figure; 
    imagesc(peakClusterLag)
    xlim([0.5 size(estX, 1)+0.5])
    ylim([0.5 size(estX, 1)+0.5])
    set(gca, 'Ytick', 1:size(estX, 1))
    
    figure;
    hold on
    for nFactor = 1:size(estX, 1)
        plot((1:1200)/4, estX(nFactor,:) + nFactor * 5);
    end
    hold off
    grid on
    
%     figure;
%     hold on
%     for nFactor = 1:size(estX, 1)
%         plot((1:1200)/4, estX(nFactor,:));
%     end
%     hold off

end

function estX         = estFactor(L, Ph, Y) % estimation of factor activity
    estX              = L'/(L*L'+diag(Ph)) * Y';
end