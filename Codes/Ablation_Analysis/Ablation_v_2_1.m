%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.1  evaluate change of number of active neurons and level of sync
%
% version using only all active neurons, independently drawn from before &
% after
% synchronization level calculated only on active neurons 
% (0 for before, activeThres for after)
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
% 

addpath('../Func');
setDir;

pThres = 0.05;
nFileList = 25:2:58;

fishListCutA = [1, 2, 3, 4, 7]; % anterior cut
fishListCutM = [8, 9, 10, 12]; % middle cut
fishListCutP = [13, 15, 16, 17]; % posterior cut

fishList = [fishListCutA, fishListCutM, fishListCutP]; % anterior cut
correlationAP = zeros(numel(fishList), 2, 2);
fracActNeuron = zeros(numel(fishList), 4, 2);
avgCorr = zeros(numel(fishList), 4, 2);
numCells = zeros(numel(fishList), 4, 2);
numActCells = zeros(numel(fishList), 4, 2);
percOverlap = zeros(numel(fishList), 1);

tagExt = '_AllActive';


for i = 1:numel(fishList)
    nFish = fishList(i);
    for nExp = 1:2
        nFile = nFileList(nFish) + nExp - 1;
        fileDirName   = fileDirNames{nFile}; %#ok<*USENS>
        fileName      = fileNames{nFile};
        
        dirImageData  = [fileDirName '/'];
        load([dirImageData, 'profile.mat'], 'segAblation');
        
        if nExp == 1
            load([tempDatDir, fileName, '.mat'], 'dff', 'activeNeuronMat', 'new_x', 'new_y');
            activeTag = activeNeuronMat>0;
        else
            load([tempDatDir, fileName, '_full.mat'], 'dff', 'activeNeuronMat', 'activeThres', 'new_x', 'new_y');
            activeTag = sum(activeNeuronMat, 2)>activeThres*10;
        end
        ind{1} = new_x<segAblation(1) & new_y<0;
        ind{2} = new_x<segAblation(1) & new_y>0;
        ind{3} = new_x>segAblation(2) & new_y<0;
        ind{4} = new_x>segAblation(2) & new_y>0;
        
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
plot([1, 2], correlationAPCombined', 'ok');
plot([1, 2], correlationAPCombined', 'k');
xlim([0.5, 2.5]);
scatter([1, 2], mean(correlationAPCombined, 1), 'or', 'linew', 3);
plot([1, 2], mean(correlationAPCombined, 1), 'r', 'linew', 3);
[h, p] = ttest(correlationAPCombined(: , 1), correlationAPCombined(: , 2));
ylim([-0.2, 1]);
set(gca, 'Xtick', [1, 2], 'Xticklabel', {'before', 'after'});
if h
    text(1.5, max(correlationAPCombined(:))*1.2, '*');
end
hold off
title('correlation between A and P')
setPrint(8, 6, [plotDir '/CorrelationAP' tagExt], 'pdf');


fishListType = {1:numel(fishListCutA), numel(fishListCutA)+1:numel(fishListCutA)+numel(fishListCutM), numel(fishListCutA)+numel(fishListCutM)+1:numel(fishList); ...
    'Cut A', 'Cut M', 'Cut P'};
figure('Position', [0, 0, 400*3, 300*2]);
hFracActAblation   = nan(size(fishListType, 2), 1);
pFracActAblaiton   = nan(size(fishListType, 2), 1);
hSyncLevelAblation    = nan(size(fishListType, 2), 1);
pSyncLevelAblation    = nan(size(fishListType, 2), 1);
for nExpType = 1:size(fishListType, 2)
    currfishList = fishListType{1, nExpType};
    fracActNeuronCombined = [fracActNeuron(currfishList, [1, 3], :); fracActNeuron(currfishList, [2, 4], :)];
    avgCorrCombined = [avgCorr(currfishList, [1, 3], :); avgCorr(currfishList, [2, 4], :)];
    subplot(2, size(fishListType, 2), nExpType);
    hold on,
    for nSec = 1:2
        xSeq = repmat([2*nSec-1; 2*nSec], 1, numel(currfishList)*2);
        ySeq = squeeze(fracActNeuronCombined(:, nSec, :))';
        scatter(xSeq(:), ySeq(:), 'ok');
        plot(xSeq, ySeq, 'k');
        [h, p] = ttest(ySeq(1, :), ySeq(2, :));
        if p < pThres
            text(2*nSec-0.5, max(ySeq(:))*1.2, '*');
        end
        scatter([2*nSec-1, 2*nSec], nanmean(squeeze(fracActNeuronCombined(:, nSec, :)))', 'or', 'linew', 3);
        plot([2*nSec-1, 2*nSec], nanmean(squeeze(fracActNeuronCombined(:, nSec, :)))', 'r', 'linew', 3);
    end
    hold off
    xlim([0, 5]);
    ylim([0, 0.5]);
    set(gca, 'Xtick', 1.5:2:5.5, 'XtickLabel', {'A', 'P'});
    title(['fraction of active neurons - ' fishListType{2, nExpType}])
    fracActRatio = squeeze(fracActNeuronCombined(:, :, 2)./fracActNeuronCombined(:, :, 1));
    [hFracActAblation(nExpType), pFracActAblaiton(nExpType)] = ttest(fracActRatio(:, 1), fracActRatio(:, 2));
    
    subplot(2, size(fishListType, 2), size(fishListType, 2)+nExpType);
    hold on,
    for nSec = 1:2
        xSeq = repmat([2*nSec-1; 2*nSec], 1, numel(currfishList)*2);
        ySeq = squeeze(avgCorrCombined(:, nSec, :))';
        scatter(xSeq(:), ySeq(:), 'ok');
        plot(xSeq, ySeq, 'k');
        [h, p] = ttest(ySeq(1, :), ySeq(2, :));
        if p < pThres
            text(2*nSec-0.5, max(ySeq(:))*1.2, '*');
        end
        scatter([2*nSec-1, 2*nSec], nanmean(squeeze(avgCorrCombined(:, nSec, :)))', 'or', 'linew', 3);
        plot([2*nSec-1, 2*nSec], nanmean(squeeze(avgCorrCombined(:, nSec, :)))', 'r', 'linew', 3);
    end
    hold off
    xlim([0, 5]);
    ylim([-0.1, 1])
    set(gca, 'Xtick', 1.5:2:5.5, 'XtickLabel', {'A', 'P'});
    title(['level of synchronization - '  fishListType{2, nExpType}])
    syncLevelRatio = squeeze(avgCorrCombined(:, :, 2)./avgCorrCombined(:, :, 1));
    [hSyncLevelAblation(nExpType), pSyncLevelAblation(nExpType)] = ttest(syncLevelRatio(:, 1), syncLevelRatio(:, 2));
end


setPrint(8*3, 6*2, [plotDir '/AblationTypeSummary' tagExt], 'pdf');

close all
