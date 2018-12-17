%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2X plot the evolution of factor center
%
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
addpath('../Func')
setDir;

nFile = 3;
exampleTimeWindow =  [30, 56, 85, 98, 122, 148];


fileName = fileNames{nFile};

load([tempDatDir, fileName, '.mat'], 'dff', 'sideSplitter', 'side', 'tracks', 'timePoints', 'activeNeuronMat', 'new_x', 'new_y', 'new_z', 'timeStep');
load([tempDatDir, 'LONOLoading_' fileName, '.mat'], 'CorrectedLMat')

if ~exist('new_x', 'var'); return; end

x                 = new_x;
y                 = new_y;
z                 = new_z;
numTime           = length(CorrectedLMat);
numNeuron         = length(side);

[~, neworder]     = sort(x);
neworder = [neworder(side(neworder)==1); neworder(side(neworder)==2)];

mColor = cbrewer('qual', 'Dark2',  8, 'cubic');
mColor            = [mColor; cbrewer('qual', 'Set2',  128, 'cubic')];
preLMat           = nan(numNeuron, 1);

linew = 1.25;
ms = 9;

z = z/max(z) * 1.8;
y = y/2;
f = figure('units', 'pixels', 'outerposition', [0 0 , 1500, 500]); hold on
set(gcf, 'PaperPositionMode', 'Auto');
whitebg('w');
set(gcf, 'Color', 'w');
set(gcf, 'InvertHardcopy', 'off');
% figure(1), hold on
% figure(2), hold on
for period = 1:numel(timePoints)
    timeRange = timePoints(period)+1:timePoints(period)+timeStep;
    radius = 0.2;
    
    LMat          = CorrectedLMat{period};
    LMat(isnan(LMat)) = 0;
    LMat(:, sum(LMat, 1)==0) = [];
    activeTag = activeNeuronMat(:, period);
    
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
            end
        end
        for nFactor = 1:size(LMat, 2)
            neuronFactor = LMat(:, nFactor)>0;
            dominateSide = ceil(mean(side(neuronFactor)) - 0.5);
            otherSide    = 3 - dominateSide;
            dominateSideNeuron = neuronFactor & side == dominateSide;
            otherSideNeuron    = neuronFactor & side == otherSide;
            xCenter = mean(x(dominateSideNeuron));
            yCenter = mean(y(dominateSideNeuron));
%             scatter(timeRange(1)/(4*3600)+17.5, xCenter, sum(dominateSideNeuron)*10, mColor(factorIndex(nFactor), :), 'filled');
            scatter3(timeRange(1)/(4*3600)+17.5, xCenter, yCenter, sum(dominateSideNeuron)*3, mColor(factorIndex(nFactor), :), 'filled');
        end
    end
%     if ismember(period, exampleTimeWindow)
%         plot([timeRange(1), timeRange(1)]/(4*3600)+17.5, [0, ceil(max(x))], 'k');
%         figure('units', 'pixels', 'outerposition', [0 0 , 300, 900]);
%         set(gcf, 'PaperPositionMode', 'Auto');
%         whitebg('w');
%         set(gcf, 'Color', 'w');
%         set(gcf, 'InvertHardcopy', 'off');
%         % dorsal view
%         hold on
%         if size(LMat,2) >= 1
%             for nFactor = 1:size(LMat, 2)
%                 neuronFactor = LMat(:, nFactor)>0;
%                 if sum(neuronFactor) <=1
%                     continue;
%                 end
%                 if length(unique(y(neuronFactor)>0))==1
%                     CHPoints = smoothedBoundary(y(neuronFactor), x(neuronFactor), radius);
%                     patch(CHPoints(:,1), CHPoints(:,2), mColor(factorIndex(nFactor), :), 'facealpha', 0.8, 'edgecolor', 'none');
%                 else
%                     if sum(y(neuronFactor)>0) == 1 || sum(y(neuronFactor)<0) == 1
%                         plot(y(neuronFactor), x(neuronFactor), 'o', 'Color', mColor(factorIndex(nFactor), :), 'MarkerFaceColor', mColor(factorIndex(nFactor), :), 'linewidth', linew, 'MarkerSize', ms)
%                     else
%                         plot(y(neuronFactor), x(neuronFactor), 's', 'Color', mColor(factorIndex(nFactor), :), 'MarkerFaceColor', mColor(factorIndex(nFactor), :), 'MarkerSize', 10, 'linewidth', linew)
%                     end
%                 end
%             end
%         end
%         plot(y(~activeTag), x(~activeTag), 'ok', 'linewidth', linew, 'MarkerFaceColor', 'w', 'MarkerSize', ms);
%         plot(y(activeTag), x(activeTag), 'ok', 'MarkerFaceColor','k', 'linewidth', linew, 'MarkerSize', ms);
%         
%         for i = 1:7
%             plot([-2, 2], [i, i], '--k');
%         end
%         ylim([0 8]);
%         xlim([-1 1]);
%         plot([0 0], [0 8], '--k');
%         set(gca, 'YDir', 'reverse');
%         set(gca, 'xticklabel', '');
%         set(gca, 'xtick', []);
%         set(gca, 'yticklabel', '');
%         set(gca, 'ytick', []);
%         box on
%         print([plotDir 'Figure_2X FA Window_' num2str(period) '.pdf'], '-dpdf', '-r0');
%         close
%     end
end
% figure(1)
set(gca, 'ydir', 'reverse');
set(gca, 'xTick', 18:20);
xlim([17.5, 21]);
set(gca, 'yTick', 0:ceil(max(x)));
set(gca, 'TickDir', 'out')
% print([plotDir 'Figure_2X FA AP Evol.pdf'], '-dpdf', '-r0');
