%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.1  (double cut) evaluate change of patterned activity based on pairwise correlation
%
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%

addpath('../Func');
setDir;

pThres = 0.05;
nFileList = 25:2:76;

fishList = [11,19,20,21,22,23,24,25,26];

numCells = zeros(numel(fishList), 6, 2);
numActCells = zeros(numel(fishList), 6, 2);
avgCorr = zeros(numel(fishList), 6, 2);
fracActNeuron = zeros(numel(fishList), 6, 2);
percOverlap = zeros(numel(fishList), 1);
correlationAP = zeros(numel(fishList), 2, 2);
correlationAM = zeros(numel(fishList), 2, 2);
correlationPM = zeros(numel(fishList), 2, 2);

tagExt = '_Double_Corr';


for i = 1:numel(fishList)
    nFish = fishList(i);
    for nExp = 1:2
        nFile = nFileList(nFish) + nExp - 1;
        fileDirName   = fileDirNames{nFile}; %#ok<*USENS>
        fileName      = fileNames{nFile};
        
        dirImageData  = [fileDirName '/'];
        load([dirImageData, 'profile.mat'], 'segAblation_A', 'segAblation_P');
        load([tempDatDir, fileName, '.mat'], 'dff', 'activeNeuronMat', 'new_x', 'new_y');

        
        if nExp == 1
            activeTag = activeNeuronMat>0;
        else
            activeTag = sum(activeNeuronMat, 2)>size(activeNeuronMat, 2)/2; %activeThres*10;
        end
        ind{1} = new_x<segAblation_A(1) & new_y<0;
        ind{2} = new_x<segAblation_A(1) & new_y>0;
        ind{3} = new_x>segAblation_A(2) & new_x<segAblation_P(1) & new_y<0;
        ind{4} = new_x>segAblation_A(2) & new_x<segAblation_P(1) & new_y>0;
        ind{5} = new_x>segAblation_P(2) & new_y<0;
        ind{6} = new_x>segAblation_P(2) & new_y>0;
        
        corMatrix = corr(dff');
        % 1. Correlaiton between average traces A and P to cutting site
        for nSec = 1:2
            % method 1.1: use average of pairwise correlation
            corMatSec = corMatrix(activeTag & ind{nSec}, activeTag & ind{nSec+4});
            correlationAP(i, nSec, nExp) = nanmean(corMatSec(:));
            corMatSec = corMatrix(activeTag & ind{nSec}, activeTag & ind{nSec+2});
            correlationAM(i, nSec, nExp) = nanmean(corMatSec(:));
            corMatSec = corMatrix(activeTag & ind{nSec+2}, activeTag & ind{nSec+4});
            correlationPM(i, nSec, nExp) = nanmean(corMatSec(:));
            % method 1.2: use correlation of average trace
            %             correlationAP(i, nSec, nExp) = corr(mean(dff(activeTag & ind{nSec}, :), 1)', mean(dff(activeTag & ind{nSec+2}, :), 1)');
        end
        % 2. fraction of active neruons in the 4 sectors
        for nSec = 1:6
            numActCells(i, nSec, nExp) = sum(activeTag & ind{nSec});
            numCells(i, nSec, nExp) = sum(ind{nSec});
            fracActNeuron(i, nSec, nExp) = sum(activeTag & ind{nSec})/sum(ind{nSec});
        end
        % 3. Average correlation between active neruons in the 4 sectors
        for nSec = 1:6
            corMatSec = corMatrix(activeTag & ind{nSec}, activeTag & ind{nSec});
            colvect = corMatSec(~tril(ones(size(corMatSec))));
            avgCorr(i, nSec, nExp) = nanmean(colvect);
        end
    end
end
% invalidExp = nan(size(avgCorr, 1), size(avgCorr, 2));
% invalidExp(~isnan(avgCorr(:, :, 1)) & ~isnan(avgCorr(:, :, 2))) = 1;
% avgCorr(:, :, 1) = avgCorr(:, :, 1) .* invalidExp;
% avgCorr(:, :, 2) = avgCorr(:, :, 2) .* invalidExp;
correlationSec = cell(3, 1);
correlationSec{1} = squeeze(nanmean(correlationAP, 2));
correlationSec{2} = squeeze(nanmean(correlationAM, 2));
correlationSec{3} = squeeze(nanmean(correlationPM, 2));



figure,
hold on
for i = 1:3
    corrCombined = correlationSec{i};
    plot([2*i-1, 2*i], corrCombined', 'ok');
    plot([2*i-1, 2*i], corrCombined', 'k');
    scatter([2*i-1, 2*i], nanmean(corrCombined, 1), 'or', 'linew', 3);
    plot([2*i-1, 2*i], nanmean(corrCombined, 1), 'r', 'linew', 3);
    p = signrank(corrCombined(: , 1), corrCombined(: , 2));
    if p < pThres
        text(2*i-0.5, max(corrCombined(:))*1.2, '*');
    end
end
ylim([-0.2, 1]);
xlim([0.5, 6.5]);

set(gca, 'Xtick', 1.5:2:5.5, 'Xticklabel', {'A-P', 'A-M', 'P-M'});
hold off
title('correlation between sectors')
setPrint(8, 6, [plotDir '/CorrelationAP' tagExt], 'pdf');


figure('Position', [0, 0, 400*2, 300]);
fracActNeuronCombined = [fracActNeuron(:, [1, 3, 5], :); fracActNeuron(:, [2, 4, 6], :)];
avgCorrCombined = [avgCorr(:, [1, 3, 5], :); avgCorr(:, [2, 4, 6], :)];
subplot(1, 2, 1);
hold on,
for nSec = 1:3
    xSeq = repmat([2*nSec-1; 2*nSec], 1, numel(fishList)*2);
    ySeq = squeeze(fracActNeuronCombined(:, nSec, :))';
    scatter(xSeq(:), ySeq(:), 'ok');
    plot(xSeq, ySeq, 'k');
    p = signrank(ySeq(1, :), ySeq(2, :));
    if p < pThres
        text(2*nSec-0.5, max(ySeq(:))*1.2, '*');
    end
    scatter([2*nSec-1, 2*nSec], nanmean(squeeze(fracActNeuronCombined(:, nSec, :)))', 'or', 'linew', 3);
    plot([2*nSec-1, 2*nSec], nanmean(squeeze(fracActNeuronCombined(:, nSec, :)))', 'r', 'linew', 3);
end
hold off
xlim([0, 7]);
ylim([-0.2, 0.8]);
set(gca, 'Xtick', 1.5:2:5.5, 'XtickLabel', {'A', 'M', 'P'});
ylabel('fraction of active neurons')
fracActRatio = squeeze(fracActNeuronCombined(:, :, 2)./fracActNeuronCombined(:, :, 1));
pFracActAblation.AP = signrank(fracActRatio(:, 1), fracActRatio(:, 3));
pFracActAblation.AM = signrank(fracActRatio(:, 1), fracActRatio(:, 2));
pFracActAblation.MP = signrank(fracActRatio(:, 2), fracActRatio(:, 3));

subplot(1, 2, 2);
hold on,
for nSec = 1:3
    xSeq = repmat([2*nSec-1; 2*nSec], 1, numel(fishList)*2);
    ySeq = squeeze(avgCorrCombined(:, nSec, :))';
    scatter(xSeq(:), ySeq(:), 'ok');
    plot(xSeq, ySeq, 'k');
    [h, p] = signrank(ySeq(1, :), ySeq(2, :));
    if p < pThres
        text(2*nSec-0.5, max(ySeq(:))*1.2, '*');
    end
    scatter([2*nSec-1, 2*nSec], nanmean(squeeze(avgCorrCombined(:, nSec, :)))', 'or', 'linew', 3);
    plot([2*nSec-1, 2*nSec], nanmean(squeeze(avgCorrCombined(:, nSec, :)))', 'r', 'linew', 3);
end
hold off
xlim([0, 7]);
ylim([-0.1, 1.2])
set(gca, 'Xtick', 1.5:2:5.5, 'XtickLabel', {'A', 'M', 'P'});
ylabel('level of synchronization')
syncLevelRatio = squeeze(avgCorrCombined(:, :, 2)./avgCorrCombined(:, :, 1));
pSyncLevelAblation.AP = signrank(syncLevelRatio(:, 1), syncLevelRatio(:, 3));
pSyncLevelAblation.AM = signrank(syncLevelRatio(:, 1), syncLevelRatio(:, 2));
pSyncLevelAblation.MP = signrank(syncLevelRatio(:, 2), syncLevelRatio(:, 3));
setPrint(8*3, 6*2, [plotDir '/AblationTypeSummary' tagExt], 'pdf');

close all