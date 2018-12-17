%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Movie S8: factor analysis result of all fish (dorsal view), dark
% background
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org



datasets = [3, 4, 7, 12, 10, 15, 16]; % not sure to include 13 or not
% datasets = [12];
noLabel = 1;
addpath('../Func');
setDir;

timeStats = zeros(numel(datasets), 3); % start, end, offset
for nFish = 1:numel(datasets)
    nFile = datasets(nFish);
    fileName = fileNames{nFile};
    load([tempDatDir, fileName, '.mat'], 'timePoints');
    load([tempDatDir, 'Leader_' fileName, '.mat'], 'timeOffset');
    timeStats(nFish, :) = [1, numel(timePoints), round(timeOffset*60)];
end
timeStats(:, 1:2) = timeStats(:, 1:2) - timeStats(:, [3, 3]);
tOffsets = timeStats(:, 1) - min(timeStats(:, 1)); % start plotting from this figure

frameW = 300;
frameH = 900;
videoStack = zeros(frameH, frameW*numel(datasets), max(timeStats(:, 2))- min(timeStats(:, 1))+1,3,'uint8')*240;
radius = 0.2;

fig = figure('units', 'pixels', 'position', [0 0 , frameW, frameH]);
set(0, 'defaultaxeslayer', 'top')
whitebg('black');
set(gcf, 'Color', 'black');
set(0 ,'defaultfigurecolor', 'black')

for nFish = 1:numel(datasets)
    
    nFile = datasets(nFish);
    fileName          = fileNames{nFile};
    load([tempDatDir, fileName, '.mat'], 'sideSplitter', 'side', 'tracks', 'timePoints', 'activeNeuronMat', 'new_x', 'new_y', 'timeStep');
    load([tempDatDir, 'LONOLoading_' fileName, '.mat'], 'CorrectedLMat')
    disp(['Processing dataset #' num2str(nFile)]);
    
    
    x                 = new_x;
    y                 = new_y;
    numNeuron         = length(side);
    
    
    mColor = cbrewer('qual', 'Dark2',  8, 'cubic');
    mColor = mColor([1 3 2 4:8], :);
    mColor            = [mColor; cbrewer('qual', 'Set3',  50, 'cubic')];
    
    
    preLMat           = nan(numNeuron, 1);
    
    
    linew = 1.25;
    y = y/2;
    
    factorIndexAll = zeros(size(mColor, 1), numel(timePoints));
    LMatAll = cell(numel(timePoints), 1);
    
    for period = 1:numel(timePoints)
        LMat          = CorrectedLMat{period};
        LMat(isnan(LMat)) = 0;
        LMat(:, sum(LMat, 1)==0) = [];
        
        % code to drop overlapped factors
        LMatNeuron    = LMat>0;
        numFactor     = size(LMatNeuron, 2);
        for nFactor   = 1:size(LMatNeuron, 2)
            if sum(LMatNeuron(:,nFactor)) >0 % skip if nFactor is dropped
                for mFactor = nFactor+1:size(LMatNeuron, 2)
                    if sum(LMatNeuron(:,mFactor)) >0 % skip if mFactor is dropped
                        if all(ismember(find(LMatNeuron(:,nFactor))', find(LMatNeuron(:,mFactor))')) % if nFactor inside any other factors, drop nFactor
                            LMatNeuron(:,nFactor) = 0;
                            continue;
                        elseif all(ismember(find(LMatNeuron(:,mFactor))', find(LMatNeuron(:,nFactor))')) % if mFactor inside any other factors, drop nFactor
                            LMatNeuron(:,mFactor) = 0;
                            continue;
                        end
                    end
                end
            end
        end
        LMat(LMatNeuron == 0)    = 0;
        LMat(:, sum(LMat, 1)==0) = []; % drop factors with zero weight
        % end of drop overlapped factor code
        
        % determine the factor index
        % break double sided factors
        LMatLeft                 = zeros(size(LMat));
        LMatRight                = zeros(size(LMat));
        LMatLeft(side == 1, :)   = LMat(side == 1, :);
        LMatRight(side == 2, :)  = LMat(side == 2, :);
        LMat                     = [LMatLeft, LMatRight];
        
        % remove the single-unit factor in movie
        LMat(:, sum(LMat>0, 1)<=1) = []; % drop factors with zero weight
        
        if size(LMat,2) >= 1
            if sum(~isnan(preLMat(:))) == 0
                factorIndex  = 1:size(LMat, 2);
                preLMat      = LMat;
%                 factorIndexAll(factorIndex, period) = factorIndex;
            else
                
                [~, maxFactorPerNeuronIndex] = max(LMat(sum(LMat, 2)>0, :), [], 2);
                sideRemoveList  = histc(maxFactorPerNeuronIndex, 1:size(LMat, 2)) <2; % remove the factor has no dominate factors
                
                LMat(:, sideRemoveList) = [];
                LMat            = LMat > 0;
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
%                     factorIndexAll(factorIndex(nFactor), period) = nFactor;
                end
            end
        end
        for n = 1:numel(factorIndex)
            factorIndexAll(factorIndex(n), period) = n;
        end
        LMatAll{period} = LMat;
    end
    % sort factor index by factor lifetime
    factorIndexAll = factorIndexAll(sum(factorIndexAll, 2)>0, :);
    [~, cid] = sort(sum(factorIndexAll>0, 2), 'descend');
    finalIndex = factorIndexAll(factorIndexAll(:, end)>0, end);
    finalColorIndex = find(factorIndexAll(:, end)>0);
    
    if mode(side(LMat(:, finalIndex(finalColorIndex==cid(1))))) == 2
        cid([1, 2]) = cid([2, 1]);
    end
    factorIndexAll = factorIndexAll(cid, :);
    
    for period = 1:numel(timePoints)
        factorIndex = factorIndexAll(factorIndexAll(:, period)>0, period);
        colorIndex  = find(factorIndexAll(:, period)>0);
        LMat = LMatAll{period};
        activeTag = activeNeuronMat(:, period);
        
        % dorsal view
        clf reset
        hold on
        plot(y(activeTag), x(activeTag), 'ow', 'MarkerFaceColor','w', 'linewidth', linew);
        if size(LMat,2) >= 1
            LMat = LMat(:, factorIndex);
            for nFactor = 1:size(LMat, 2)
                neuronFactor = LMat(:, nFactor)>0;
                dominateSide = ceil(mean(side(neuronFactor)) - 0.5);
                otherSide    = 3 - dominateSide;
                dominateSideNeuron = neuronFactor & side == dominateSide;
                otherSideNeuron    = neuronFactor & side == otherSide;
                CHPoints = smoothedBoundary(y(dominateSideNeuron), x(dominateSideNeuron), radius);
                patch(CHPoints(:,1), CHPoints(:,2), mColor(colorIndex(nFactor), :), 'facealpha', 0.6, 'edgecolor', 'none');
                if sum(otherSideNeuron)>0
                    if sum(otherSideNeuron) > 0
                        CHPoints = smoothedBoundary(y(otherSideNeuron), x(otherSideNeuron), radius);
                        patch(CHPoints(:,1), CHPoints(:,2), mColor(colorIndex(nFactor), :), 'facealpha', 0.6, 'edgecolor', 'none');
                    end
                end
            end
        end
        plot(y(~activeTag), x(~activeTag), 'ow', 'linewidth', linew, 'MarkerFaceColor', 'k');
        ylim([0 ceil(max(x))]);
        xlim([-1 1]);
        set(gca, 'YDir', 'reverse');
        set(gca, 'xticklabel', '');
        set(gca, 'xtick', []);
        set(gca, 'yticklabel', '');
        set(gca, 'ytick', []);
        if ~noLabel
            for i = 1:ceil(max(x)-1)
                plot([-2, 2], [i, i],'--w');
            end
            plot([0, 0], [0, ceil(max(x))], '--w');
            box on
        else
            axis off
        end
        
        hold off
        videoStack(1:frameH, frameW*(nFish-1) + (1:frameW), period+tOffsets(nFish),:) = imresize(frame2im(getframe(fig)), [frameH, frameW]);
    end
    for tend = period+tOffsets(nFish)+1 : size(videoStack, 3)
        videoStack(1:frameH, frameW*(nFish-1) + (1:frameW), tend,:) = imresize(frame2im(getframe(fig)), [frameH, frameW]);
    end
end

tagExt = '';
if noLabel
    tagExt = '_nolabel';
end
% frame size in pixels
writeImage(videoStack(:, :, :, 1), [plotDir '/Movie_S8' tagExt '_r.klb']);
writeImage(videoStack(:, :, :, 2), [plotDir '/Movie_S8' tagExt '_g.klb']);
writeImage(videoStack(:, :, :, 3), [plotDir '/Movie_S8' tagExt '_b.klb']);

