%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Movie S6 Covariance analysis: Factor analysis -- loadin matrix evolution --
% reordering
%
% Dark background
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org


nFile = 3;
noLabel = 0;

addpath('../Func');
setDir;
fileName          = fileNames{nFile}; 
load([tempDatDir, fileName, '.mat'], 'dff', 'sideSplitter', 'side', 'tracks', 'timePoints', 'activeNeuronMat', 'new_x', 'new_y', 'new_z', 'timeStep');
load([tempDatDir, 'LONOLoading_' fileName, '.mat'], 'CorrectedLMat')


x                 = new_x;
y                 = new_y;
z                 = new_z;
numTime           = length(CorrectedLMat);
numNeuron         = length(side);

[~, neworder]     = sort(x);
neworder = [neworder(side(neworder)==2); neworder(side(neworder)==1)];

% mColor = cbrewer('qual', 'Dark2',  8, 'cubic');
% mColor            = [mColor; cbrewer('qual', 'Set2',  128, 'cubic')];

mColor = cbrewer('qual', 'Dark2',  8, 'cubic');
mColor = mColor([1 3 2 4:8], :);
mColor            = [mColor; cbrewer('qual', 'Set3',  15, 'cubic')];

preLMat           = nan(numNeuron, 1);

if noLabel
    video          = VideoWriter([plotDir '\Movie_S6_nolabel.avi'], 'Uncompressed AVI');
else
    video          = VideoWriter([plotDir '\Movie_S6.avi'], 'Uncompressed AVI');
end
video.FrameRate = 10;
open(video);
allID = [];

% frame size in pixels
frameW = 2000;
frameH = 1000;
fig = figure('units', 'pixels', 'position', [0 0 , frameW, frameH]);
set(0, 'defaultaxeslayer', 'top')
whitebg('black');
set(gcf, 'Color', 'black');
set(0 ,'defaultfigurecolor', 'black')

linew = 1.25;
z = z/max(z) * 1.8;
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
            factorIndexAll(factorIndex, period) = factorIndex;
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
                factorIndexAll(factorIndex(nFactor), period) = nFactor;
            end
        end
    end
    LMatAll{period} = LMat;
end

% sort factor index by factor lifetime
[~, cid] = sort(sum(factorIndexAll>0, 2), 'descend');
factorIndexAll = factorIndexAll(cid, :);

for period = 1:numel(timePoints)
    factorIndex = factorIndexAll(factorIndexAll(:, period)>0, period);
    colorIndex  = find(factorIndexAll(:, period)>0);
    LMat = LMatAll{period};
    activeTag = activeNeuronMat(:, period);
    LMat = LMat(:, factorIndex);
    
    timeRange = timePoints(period)+1:timePoints(period)+timeStep;
    clf reset
    radius = 0.2;    
    
    % plot calcium traces
    subplot(3, 2, [1 3 5]);
    hold on
    for i = 1:size(LMat, 1)
        if ~activeTag(i)
            plot(linspace(0, 5, numel(timeRange)), zscore(dff(i, timeRange))+find(neworder==i)*4, 'Color', [.2, .2, .2], 'linewidth', linew);
        end
    end
    for i = 1:size(LMat, 1)
        if activeTag(i)
            if ~any(isnan(LMat(i, :))) && sum(LMat(i, :))>0 && size(LMat,2) >= 1
                [~, nFactor] = max(LMat(i, :));
                plot(linspace(0, 5, numel(timeRange)), zscore(dff(i, timeRange))+find(neworder==i)*4, 'Color', mColor(colorIndex(nFactor), :), 'linewidth', linew);
            else
                plot(linspace(0, 5, numel(timeRange)), zscore(dff(i, timeRange))+find(neworder==i)*4, 'w', 'linewidth', linew);
            end
        end
    end
    xlim([0, 5])
    ylim([0, size(LMat, 1)*4+10]);
    set(gca, 'YTickLabel', '');
    if ~noLabel
        plot([0, 5], [sum(side==1) * 4, sum(side==1) * 4], 'w--');
    else
        axis off
    end
    hold off
    
    % plot spatial organization
    % dorsal view
    subplot(3, 2, 2)
    hold on
    plot(x(activeTag), y(activeTag), 'ow', 'MarkerFaceColor','w', 'linewidth', linew);
    if size(LMat,2) >= 1
        for nFactor = 1:size(LMat, 2)
            neuronFactor = LMat(:, nFactor)>0;
            dominateSide = ceil(mean(side(neuronFactor)) - 0.5);
            otherSide    = 3 - dominateSide;
            dominateSideNeuron = neuronFactor & side == dominateSide;
            otherSideNeuron    = neuronFactor & side == otherSide;
            CHPoints = smoothedBoundary(x(dominateSideNeuron), y(dominateSideNeuron), radius);
            patch(CHPoints(:,1), CHPoints(:,2), mColor(colorIndex(nFactor), :), 'facealpha', 0.6, 'edgecolor', 'none');
            if sum(otherSideNeuron)>0
                if sum(otherSideNeuron) > 0
                    CHPoints = smoothedBoundary(x(otherSideNeuron), y(otherSideNeuron), radius);
                    patch(CHPoints(:,1), CHPoints(:,2), mColor(colorIndex(nFactor), :), 'facealpha', 0.6, 'edgecolor', 'none');
                end
            end
        end
    end
    plot(x(~activeTag), y(~activeTag), 'ow', 'linewidth', linew, 'MarkerFaceColor', 'k');
    if ~noLabel
        for i = 1:ceil(max(x))
            plot([i, i], [-2, 2], '--w');
        end
    else
        axis off
    end
    xlim([0 ceil(max(x))]);
    ylim([-1 1]);
    %         text(ceil(max(x))+0.5, 0.4, num2str(period), 'FontSize', 24)
    set(gca, 'YTicklabel', '');

    hold off
    
    % side view - left
    subplot(3, 2, 6)
    hold on
    plot(x(activeTag & side==1), z(activeTag & side==1), 'ow', 'MarkerFaceColor','w', 'linewidth', linew);
    if size(LMat,2) >= 1
        for nFactor = 1:size(LMat, 2)
            neuronFactor = LMat(:, nFactor)>0;
            if sum(neuronFactor & side==1)>0
                CHPoints = smoothedBoundary(x(neuronFactor & side==1), z(neuronFactor & side==1), radius);
                patch(CHPoints(:,1), CHPoints(:,2), mColor(colorIndex(nFactor), :), 'facealpha', 0.6, 'edgecolor', 'none');
            end
        end
    end
    plot(x(~activeTag & side==1), z(~activeTag & side==1), 'ow', 'linewidth', linew, 'MarkerFaceColor', 'k');
    if ~noLabel
        for i = 1:ceil(max(x))
            plot([i, i], [0, 2], '--w');
        end
    else
        axis off
    end
    xlim([0 ceil(max(x))]);
    ylim([0 2]);    
    set(gca, 'YTicklabel', '');

    hold off
    
    % side view - right
    subplot(3, 2, 4)
    hold on
    plot(x(activeTag & side==2), z(activeTag & side==2), 'ow', 'MarkerFaceColor','w', 'linewidth', linew);
    if size(LMat,2) >= 1
        for nFactor = 1:size(LMat, 2)
            neuronFactor = LMat(:, nFactor)>0;
            if sum(neuronFactor & side==2)>0
                CHPoints = smoothedBoundary(x(neuronFactor & side==2), z(neuronFactor & side==2), radius);
                patch(CHPoints(:,1), CHPoints(:,2), mColor(colorIndex(nFactor), :), 'facealpha', 0.6, 'edgecolor', 'none');
            end
        end
    end
    plot(x(~activeTag & side==2), z(~activeTag & side==2), 'ow', 'linewidth', linew, 'MarkerFaceColor', 'k');
    if ~noLabel
        for i = 1:ceil(max(x))
            plot([i, i], [0, 2], '--w');
        end
    else
        axis off
    end
    xlim([0 ceil(max(x))]);
    ylim([0 2]);
    set(gca, 'YTicklabel', '');
    hold off
    frame = getframe(fig);
    writeVideo(video, frame);
    
    allID = unique([allID; colorIndex]);
end
close(video);
close;



    

