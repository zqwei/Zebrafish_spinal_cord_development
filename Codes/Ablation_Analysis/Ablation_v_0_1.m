%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0.  evaluate change of number of active neurons and level of sync
%
% plot on dark background
% 
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
% 

addpath('../Func');
setDir;

nFileList = 25:2:58;

fishListCutA = [1, 2, 3, 4, 7]; % anterior cut
fishListCutM = [5, 8, 9, 10, 12]; % middle cut
fishListCutP = [13, 15, 16, 17]; % posterior cut


fishList = [fishListCutA, fishListCutM, fishListCutP]; % anterior cut
correlationAP = zeros(numel(fishList), 2, 2);
fracActNeuron = zeros(numel(fishList), 4, 2);
avgCorr = zeros(numel(fishList), 4, 2);
numCells = zeros(numel(fishList), 4, 2);
numActCells = zeros(numel(fishList), 4, 2);
for i = 1:numel(fishList)
    nFish = fishList(i);
    for nExp = 1:2
        nFile = nFileList(nFish) + nExp - 1;
        fileDirName   = fileDirNames{nFile}; %#ok<*USENS>
        fileName      = fileNames{nFile};
        
        dirImageData  = [fileDirName '/'];
        load([dirImageData, 'profile.mat'], 'segAblation');
        load ([tempDatDir, fileName, '.mat'], 'dff', 'activeNeuronMat', 'new_x', 'new_y', 'new_z');
        
        
        activeThres = floor(size(activeNeuronMat, 2)/4); %1
        activeTag = sum(activeNeuronMat, 2)> activeThres;
        ind{1} = new_x<segAblation(1) & new_y<0;
        ind{2} = new_x<segAblation(1) & new_y>0;
        ind{3} = new_x>segAblation(2) & new_y<0;
        ind{4} = new_x>segAblation(2) & new_y>0;
        
        %         if ismember(nFish, fishListCutP)
        %             ind = ind([3, 4, 1, 2]);
        %         end
        
        corMatrix = corr(dff');
        % 1. Correlaiton between average traces A and P to cutting site
        for nSec = 1:2
            % method 1.1: use average of pairwise correlation
            corMatSec = corMatrix(activeTag & ind{nSec}, activeTag & ind{nSec+2});
            correlationAP(i, nSec, nExp) = nanmean(corMatSec(:));
            % method 1.2: use correlation of average trace
            %             correlationAP(i, nSec, nExp) = corr(mean(dff(activeTag & ind{nSec}, :), 1)', mean(dff(activeTag & ind{nSec+2}, :), 1)');
        end
        % 2. fraction of active neruons in the 4 sectors
        for nSec = 1:4
            numActCells(i, nSec, nExp) = sum(activeTag & ind{nSec});
            numCells(i, nSec, nExp) = sum(ind{nSec});
            fracActNeuron(i, nSec, nExp) = sum(activeTag & ind{nSec})/sum(ind{nSec});
        end
        % 3. Average correlation between active neruons in the 4 sectors
        for nSec = 1:4
            corMatSec = corMatrix(activeTag & ind{nSec}, activeTag & ind{nSec});
            colvect = corMatSec(~tril(ones(size(corMatSec))));
            avgCorr(i, nSec, nExp) = nanmean(colvect);
        end
    end
end
invalidExp = nan(size(avgCorr, 1), size(avgCorr, 2));
invalidExp(~isnan(avgCorr(:, :, 2))) = 1;
avgCorr(:, :, 1) = avgCorr(:, :, 1) .* invalidExp;
correlationAPCombined = squeeze(nanmean(correlationAP, 2));


figure,
hold on
plot([1, 2], correlationAPCombined', 'ow');
plot([1, 2], correlationAPCombined', 'w');
xlim([0.5, 2.5]);
scatter([1, 2], mean(correlationAPCombined, 1), 'or', 'linew', 3);
plot([1, 2], mean(correlationAPCombined, 1), 'r', 'linew', 3);
ylim([-0.2, 1]);
set(gca, 'Xtick', [1, 2], 'Xticklabel', {'before', 'after'});
hold off
title('correlation between A and P')
whitebg('black');
set(gcf, 'Color', 'black');
set(gcf, 'InvertHardCopy', 'off', 'PaperPositionMode', 'Auto');
print([plotDir '/CorrelationAP.pdf'], '-dpdf', '-r0');

fishListType = {1:numel(fishListCutA), numel(fishListCutA)+1:numel(fishListCutA)+numel(fishListCutM), numel(fishListCutA)+numel(fishListCutM)+1:numel(fishList); ...
    'Cut A', 'Cut M', 'Cut P'};
figure('Position', [0, 0, 400*3, 300*2]);
for nExpType = 1:size(fishListType, 2)
    currfishList = fishListType{1, nExpType};
    fracActNeuronCombined = [fracActNeuron(currfishList, [1, 3], :); fracActNeuron(currfishList, [2, 4], :)];
    avgCorrCombined = [avgCorr(currfishList, [1, 3], :); avgCorr(currfishList, [2, 4], :)];
    subplot(2, size(fishListType, 2), nExpType);
    hold on,
    for nSec = 1:2
        xSeq = repmat([2*nSec-1; 2*nSec], 1, numel(currfishList)*2);
        ySeq = squeeze(fracActNeuronCombined(:, nSec, :))';
        scatter(xSeq(:), ySeq(:), 'ow');
        plot(xSeq, ySeq, 'w');
        scatter([2*nSec-1, 2*nSec], nanmean(squeeze(fracActNeuronCombined(:, nSec, :)))', 'or', 'linew', 3);
        plot([2*nSec-1, 2*nSec], nanmean(squeeze(fracActNeuronCombined(:, nSec, :)))', 'r', 'linew', 3);
    end
    hold off
    xlim([0, 5]);
    ylim([0, 0.5]);
    set(gca, 'Xtick', 1.5:2:5.5, 'XtickLabel', {'A', 'P'});
    title(['fraction of active neurons - ' fishListType{2, nExpType}])
    %     whitebg('black');
    %     set(gcf, 'Color', 'black');
    %     set(gcf, 'InvertHardCopy', 'off', 'PaperPositionMode', 'Auto');
    
    subplot(2, size(fishListType, 2), nExpType+size(fishListType, 2));
    hold on,
    for nSec = 1:2
        xSeq = repmat([2*nSec-1; 2*nSec], 1, numel(currfishList)*2);
        ySeq = squeeze(avgCorrCombined(:, nSec, :))';
        scatter(xSeq(:), ySeq(:), 'ow');
        plot(xSeq, ySeq, 'w');
        scatter([2*nSec-1, 2*nSec], nanmean(squeeze(avgCorrCombined(:, nSec, :)))', 'or', 'linew', 3);
        plot([2*nSec-1, 2*nSec], nanmean(squeeze(avgCorrCombined(:, nSec, :)))', 'r', 'linew', 3);
    end
    hold off
    xlim([0, 5]);
    ylim([-0.1, 1])
    set(gca, 'Xtick', 1.5:2:5.5, 'XtickLabel', {'A', 'P'});
    title(['level of synchronization - '  fishListType{2, nExpType}])
end
whitebg('black');
set(gcf, 'Color', 'black');
set(gcf, 'InvertHardCopy', 'off', 'PaperPositionMode', 'Auto');
print([plotDir '/AblationTypeSummary.pdf'], '-dpdf', '-r0');
