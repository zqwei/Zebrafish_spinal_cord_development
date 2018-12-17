%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2B plot snapshots of FA result on atlas
%
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
% 

addpath('../Func')
setDir;
load([tempDatDir '/' fileNames{3} '.mat'], 'dff', 'new_x', 'new_y', 'timePoints', 'activeNeuronMat');
load([tempDatDir '/LONOLoading_' fileNames{3} '.mat'], 'CorrectedLMat');


example_windows = [30, 56, 85, 98, 122, 148];

% mColor = [cbrewer('qual', 'Set1',  7, 'cubic')];

mColor = cbrewer('qual', 'Dark2',  8, 'cubic');
mColor = mColor([1, 2, 3, 5, 4, 6, 8], :);
% mColor = mColor(randperm(7), :);
nNeurons = numel(new_x);
preLMat           = nan(nNeurons, 1);

ratio3D = 6;
linew = 1.25;
ms = 9;

x = new_x;
y = new_y/2;
for period = example_windows
    
    timeRange = timePoints(period)+1:timePoints(period)+240;
    radius = 0.3;
    
    LMat          = CorrectedLMat{period};
    LMat(isnan(LMat)) = 0;
    LMat(:, sum(LMat, 1)==0) = [];
    activeTag = activeNeuronMat(:, period);
    % drop factors completely contained within another factor
    invIndex = [];
    currLMat = LMat>0;
    for ii = 1:size(LMat, 2)
        for jj = ii+1:size(LMat, 2)
            if sum(ismember(find(currLMat(:, ii)), find(currLMat(:, jj))))==sum(currLMat(:, ii))
                invIndex = [invIndex, ii];
            elseif sum(ismember(find(currLMat(:, ii)), find(currLMat(:, jj))))==sum(currLMat(:, jj))
                invIndex = [invIndex, jj];
            end
        end
    end
    LMat(:, invIndex) = [];
    % determine the factor index
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
    end
         
    if ~ismember(period, example_windows)
        continue;
    else
        figure('units', 'pixels', 'outerposition', [0 0 , 300, 900]);
        set(gcf, 'PaperPositionMode', 'Auto');  
        whitebg('w');
        set(gcf, 'Color', 'w');
        set(gcf, 'InvertHardcopy', 'off');
        % dorsal view
        hold on
        if size(LMat,2) >= 1
            for nFactor = 1:size(LMat, 2)
                neuronFactor = LMat(:, nFactor)>0;
                if sum(neuronFactor) <=1 
                    continue;
                end
                if length(unique(y(neuronFactor)>0))==1
                    CHPoints = smoothedBoundary(y(neuronFactor), x(neuronFactor), radius);
                    patch(CHPoints(:,1), CHPoints(:,2), mColor(factorIndex(nFactor), :), 'facealpha', 0.8, 'edgecolor', 'none');
                else
                    if sum(y(neuronFactor)>0) == 1 || sum(y(neuronFactor)<0) == 1
                        plot(y(neuronFactor), x(neuronFactor), 'o', 'Color', mColor(factorIndex(nFactor), :), 'MarkerFaceColor', mColor(factorIndex(nFactor), :), 'linewidth', linew, 'MarkerSize', ms)
                    else
                        plot(y(neuronFactor), x(neuronFactor), 's', 'Color', mColor(factorIndex(nFactor), :), 'MarkerFaceColor', mColor(factorIndex(nFactor), :), 'MarkerSize', 10, 'linewidth', linew)
                    end
                end
            end
        end
        plot(y(~activeTag), x(~activeTag), 'ok', 'linewidth', linew, 'MarkerFaceColor', 'w', 'MarkerSize', ms);
        plot(y(activeTag), x(activeTag), 'ok', 'MarkerFaceColor','k', 'linewidth', linew, 'MarkerSize', ms);

        for i = 1:7
            plot([-2, 2], [i, i], '--k');
        end
        ylim([0 8]);
        xlim([-1 1]);
        plot([0 0], [0 8], '--k');
        set(gca, 'YDir', 'reverse');
        set(gca, 'xticklabel', '');
        set(gca, 'xtick', []);
        set(gca, 'yticklabel', '');
        set(gca, 'ytick', []);
        box on
        print([plotDir 'Figure_2B FA Window_' num2str(period) '.pdf'], '-dpdf', '-r0');
        close;
    end
end