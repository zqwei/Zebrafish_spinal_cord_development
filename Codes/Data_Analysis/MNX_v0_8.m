%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.  Covariance analysis: Factor analysis -- loadin matrix evolution --
% reordering
% 
% Yinan's version of FACluster_v0_5
% 
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
% 
% Integrated by Ziqiang Wei


function MNX_v0_8(nFile)            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile};  %#ok<*USENS>
    load([tempDatDir, fileName, '.mat'], 'sideSplitter', 'side', 'tracks', 'timePoints', 'activeNeuronMat', 'new_x', 'new_y', 'new_z', 'mnx'); 
    load([tempDatDir, 'EvoLoading_' fileName, '_v2.mat'], 'networkMat')
    if ~exist('mnx', 'var') || ~exist('new_x', 'var')
        return;
    end
    x                 = new_x;
    y                 = new_y;
    z                 = new_z;
    numTime           = length(networkMat);
    numNeuron         = length(side);
    
    [~, neworder]     = sort(x);
    neworder = [neworder(side(neworder)==1); neworder(side(neworder)==2)];
    
    
    % video specs
    mColor = cbrewer('qual', 'Dark2',  8, 'cubic');
    mColor = [mColor; cbrewer('qual', 'Set2',  128, 'cubic')];
    delayColor = cbrewer('seq', 'YlOrRd',  64, 'cubic');
    delayColor = flipud(delayColor);
    preLMat           = nan(numNeuron, 1);
    video          = VideoWriter([plotDir '\movie_phase_' fileName '.avi'], 'Uncompressed AVI');
    video.FrameRate = 10;
    open(video);
    frameW = 2000;
    frameH = 1000;
    fig = figure('units', 'pixels', 'position', [0 0 , frameW, frameH]);
    
    linew = 2.25;
    z = z/max(z) * 1.8;
    y = y/2;
    maxDelay = 15;
    
    for nTime = 1:numTime
        clf reset
        radius = 0.2;
        
        hold on
        % basic setting of the plots
        for i = 1:ceil(max(x))
            plot([i, i], [-2, 2], '--', 'color', [0.3 0.3 0.3]);
        end
        text(ceil(max(x)) + 0.1, 0.8, [num2str(nTime) ' min'],'fontsize', 12)
        xlim([0 ceil(max(x))+1]);
        ylim([-1 1]);
        
        networkTime    = networkMat{nTime, 1}; %#ok<*USENS>
        networkCorr    = networkMat{nTime, 2};
        numUnitsFactor = [networkTime.nUnit];
        FactorIndex    = numUnitsFactor > 1;
        FactorIndex(1) = false;
        FactorMNX      = [networkTime.mnxFrac];
        
        activeTag      = true(numNeuron, 1);
        activeTag(networkTime(1).neuronIndex) = false;
        
        plot(x(activeTag & mnx), y(activeTag & mnx), 'ok', 'linewidth', linew);
        plot(x(activeTag & ~mnx), y(activeTag & ~mnx), 'or', 'linewidth', linew);
        plot(x(~activeTag & mnx), y(~activeTag & mnx), 'ok', 'linewidth', linew, 'MarkerFaceColor', 'w');
        plot(x(~activeTag & ~mnx), y(~activeTag & ~mnx), 'or', 'linewidth', linew, 'MarkerFaceColor', 'w');
        
        
        if sum(FactorIndex) == 0
            continue;
        end
        
        LMat                = false(numNeuron, sum(FactorIndex));
        FactorIndexSet      = find(FactorIndex);
        for nFactor         = 1:length(FactorIndexSet)
            mFactor         = FactorIndexSet(nFactor);
            LMat(networkTime(mFactor).neuronIndex, nFactor) = true;
        end
        
        % determine the factor index
        if sum(~isnan(preLMat(:))) == 0
            factorIndex  = 1:size(LMat, 2);
            preLMat      = LMat;
        else
%             [~, maxFactorPerNeuronIndex] = max(LMat(sum(LMat, 2)>0, :), [], 2);
%             sideRemoveList  = histc(maxFactorPerNeuronIndex, 1:size(LMat, 2)) <2; % remove the factor has no dominate factors
%             LMat(:, sideRemoveList) = [];
            sizeLMat        = sum(LMat, 1);
            [~, indexLMat]  = sort(sizeLMat, 'descend');
            LMat            = LMat(:, indexLMat);
            factorIndex     = zeros(size(LMat, 2), 1);
            % compute similarity matrix
            similarityScore = zeros(size(LMat, 2), size(preLMat, 2));
            for nFactor     = 1:size(LMat, 2)
                if sum(isnan(LMat(:)))>0; keyboard();end
                similarityScore(nFactor, :) = sum(bsxfun(@and, LMat(:, nFactor), preLMat));
            end
            % check if any prefactor has no connection with new factors
            % decide which factor is not included in preLMatIndex
            % maxIndex is the factor with the maximum coverage with the prefactors
            [~, maxIndex]   = max(similarityScore, [], 1);
            % check if any prefactor is merged (factor has maximum coverages with more than one prefactors, pick the larger one as its index)
            for nFactor     = 1:size(LMat, 2)
                nFacotrNumPreFactor = sum(maxIndex == nFactor);
                switch nFacotrNumPreFactor
                    case 0
                        preLMat = [preLMat, LMat(:, nFactor)];
                        factorIndex(nFactor) = size(preLMat, 2);
                    case 1
                        if similarityScore(nFactor, maxIndex == nFactor) == 0
                            preLMat = [preLMat, LMat(:, nFactor)];
                            factorIndex(nFactor) = size(preLMat, 2);
                        else
                            preLMatIndex         = find(maxIndex == nFactor);
                            factorIndex(nFactor) = preLMatIndex;
                            preLMat(:, preLMatIndex) = preLMat(:, preLMatIndex) | LMat(:, nFactor);
                        end
                    otherwise
                        preLMatIndex         = find(maxIndex == nFactor);
                        [~, nFactorMaxIndex] = max(similarityScore(nFactor, preLMatIndex));
                        factorIndex(nFactor) = preLMatIndex(nFactorMaxIndex);
                        preLMat(:, preLMatIndex(nFactorMaxIndex)) = preLMat(:, preLMatIndex(nFactorMaxIndex)) | LMat(:, nFactor);
                end
            end
        end

        % plot of FA area and location
        for nFactor = 1:size(LMat, 2)
            neuronFactor = LMat(:, nFactor)>0;
            if length(unique(side(neuronFactor)))==1
                CHPoints = smoothedBoundary(x(neuronFactor), y(neuronFactor), radius);
                CHPoints(end+1, :) = CHPoints(1, :);
                plot(CHPoints(:,1), CHPoints(:,2), '-', 'linewid', 2.0, 'color', mColor(factorIndex(nFactor), :));
            else
                dominateSide = ceil(mean(side(neuronFactor)) - 0.5);
                otherSide    = 3 - dominateSide;
                CHPoints = smoothedBoundary(x(neuronFactor & side==dominateSide), y(neuronFactor & side==dominateSide), radius);
                CHPoints(end+1, :) = CHPoints(1, :);
                plot(CHPoints(:,1), CHPoints(:,2), '-', 'linewid', 2.0, 'color', mColor(factorIndex(nFactor), :));
                plot(x(neuronFactor & side==otherSide), y(neuronFactor & side==otherSide), '.', 'color', mColor(factorIndex(nFactor), :));
            end
        end
        
        for nFactor     = find(FactorIndex)
            for mFactor = find(FactorIndex)
                if nFactor < mFactor
                    delayTime = networkCorr.delayMat(nFactor-1, mFactor-1);
                    delayTime = abs(delayTime);
                    delayCorr = abs(networkCorr.corrMat(nFactor-1, mFactor-1));
                    isContra  = networkCorr.IpsiIndex(nFactor-1, mFactor-1) == 0;
                    if isContra
                        FALoc1_x = mean(x(networkTime(nFactor).neuronIndex));
                        FALoc2_x = mean(x(networkTime(mFactor).neuronIndex));
                        FALoc1_y = mean(y(networkTime(nFactor).neuronIndex));
                        FALoc2_y = mean(y(networkTime(mFactor).neuronIndex)); 
                        if isnan(delayTime) || delayTime>maxDelay; delayTime = maxDelay; end
                        if isnan(delayCorr); delayCorr = 0; end
                        plot([FALoc1_x, FALoc2_x], [FALoc1_y, FALoc2_y], '-', 'color', delayColor(ceil(delayTime/maxDelay*63)+1, :), 'linewidth', delayCorr*10+0.001)
                    end
                end
            end
        end

        
        hold off
        pause(0.1);
        frame        = getframe;
        writeVideo(video, frame);
        
        
    end
    
    close(video);
    close;
    
end