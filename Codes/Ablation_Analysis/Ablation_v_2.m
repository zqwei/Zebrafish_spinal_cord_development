%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  (single cut) evaluate change of level of sync  using FA result
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
% 

function Ablation_v_2(fishListCutA, fishListCutM, fishListCutP)
addpath('../Func');
setDir;

pThres = 0.05;
nFileList = 25:2:100;

fishList = [fishListCutA, fishListCutM, fishListCutP]; % anterior cut
numCells = zeros(numel(fishList), 4, 2);
numActCells = zeros(numel(fishList), 4, 2);
fracFacNeuron = zeros(numel(fishList), 4, 2);
factorSpan = zeros(numel(fishList), 4, 2);
correlationAP = zeros(numel(fishList), 2, 2);

tagExt = '_FA';
try

for i = 1:numel(fishList)
    nFish = fishList(i);
    for nExp = 1:2
        nFile = nFileList(nFish) + nExp - 1;
        fileDirName   = fileDirNames{nFile}; %#ok<*USENS>
        fileName      = fileNames{nFile};
        
        dirImageData  = [fileDirName '/'];
        load([dirImageData, 'profile.mat'], 'segAblation');
        load([tempDatDir, 'FALONO_Average_' fileName, '.mat'], 'LMat');
        load([tempDatDir, fileName, '.mat'], 'dff', 'activeNeuronMat', 'new_x', 'new_y');

        ind{1} = new_x<segAblation(1) & new_y<0;
        ind{2} = new_x<segAblation(1) & new_y>0;
        ind{3} = new_x>segAblation(2) & new_y<0;
        ind{4} = new_x>segAblation(2) & new_y>0;
        activeTag = sum(activeNeuronMat, 2)>0;
        factorTag = sum(LMat, 2) >0;
%         factorTag = sum(double(LMat>0) .* repmat(1:size(LMat, 2), size(LMat, 1), 1), 2);
        % 1. Correlaiton between average traces A and P to cutting site
        for nSec = 1:2
            % method 1.1: use average of pairwise correlation
            corMatrix = corr(dff');
            corMatSec = corMatrix(activeTag & ind{nSec}, activeTag & ind{nSec+2});
            correlationAP(i, nSec, nExp) = nanmean(corMatSec(:));
%             % method 1.2: use correlation of average trace
%             correlationAP(i, nSec, nExp) = corr(mean(dff(activeTag & ind{nSec}, :), 1)', mean(dff(activeTag & ind{nSec+2}, :), 1)');
        end
        % 2. fraction of factored neruons in the 4 sectors
        for nSec = 1:4
            numActCells(i, nSec, nExp) = sum(activeTag & ind{nSec});
            numCells(i, nSec, nExp) = sum(ind{nSec});
            if sum(factorTag & ind{nSec})>0
                fracFacNeuron(i, nSec, nExp) = sum(factorTag & ind{nSec})/sum(activeTag & ind{nSec});
                factorSpan(i, nSec, nExp) = (max(new_x(factorTag & ind{nSec})) - min(new_x(factorTag & ind{nSec})))/(max(new_x(activeTag & ind{nSec})) - min(new_x(activeTag & ind{nSec})));
            else
                 fracFacNeuron(i, nSec, nExp) = 0;
                 factorSpan(i, nSec, nExp) = 0;
            end
        end
    end
end

correlationAPCombined = squeeze(nanmean(correlationAP, 2));
figure, 
hold on
plot([1, 2], correlationAPCombined', 'ok');
plot([1, 2], correlationAPCombined', 'k');
xlim([0.5, 2.5]);
scatter([1, 2], mean(correlationAPCombined, 1), 60, 'xr', 'linew', 3);
plot([1, 2], mean(correlationAPCombined, 1), 'r', 'linew', 3);
ylim([0, 1]);
set(gca, 'Xtick', [1, 2], 'Xticklabel', {'before', 'after'});
hold off
title('correlation between A and P')
setPrint(8, 6, [plotDir '/CorrelationAP'],'pdf');

fishListType = {1:numel(fishListCutA), numel(fishListCutA)+1:numel(fishListCutA)+numel(fishListCutM), numel(fishListCutA)+numel(fishListCutM)+1:numel(fishList); ...
    'Cut A', 'Cut M', 'Cut P'};
figure('Position', [0, 0, 400*3, 300*2]);
pFracFacNeuronAblaiton = nan(size(fishListType, 2), 1);
pFactorSpanAblation    = nan(size(fishListType, 2), 1);
ratioList              = []; % fracFacNeuronRatio-A,P, fracFactorSpanRatio-A,P, expType  
ratioStats             = zeros(size(fishListType, 2), 4);% mean-A, mean-P, std-A, std-P
spanStats             = zeros(size(fishListType, 2), 4);% mean-A, mean-P, std-A, std-P of factor span
for nExpType = 1:size(fishListType, 2)
    currfishList = fishListType{1, nExpType};
    fracFacNeuronCombined = [fracFacNeuron(currfishList, [1, 3], :); fracFacNeuron(currfishList, [2, 4], :)];
    factorSpanCombined = [factorSpan(currfishList, [1, 3], :); factorSpan(currfishList, [2, 4], :)];
    subplot(2, size(fishListType, 2), nExpType);
    hold on,
    for nSec = 1:2
        xSeq = repmat([2*nSec-1; 2*nSec], 1, numel(currfishList)*2);
        ySeq = squeeze(fracFacNeuronCombined(:, nSec, :))';
        scatter(xSeq(:), ySeq(:), 'ok');
        plot(xSeq, ySeq, 'k');
        p = signrank(ySeq(1, :), ySeq(2, :));
        if p < pThres
            text(2*nSec-0.5, max(ySeq(:))*1.2, '*');
            disp(['nExpType=' num2str(nExpType) ', nSec=' num2str(nSec) ', p=' num2str(p)]);
        end
        scatter([2*nSec-1, 2*nSec], nanmean(squeeze(fracFacNeuronCombined(:, nSec, :)))', 'or', 'linew', 3);
        plot([2*nSec-1, 2*nSec], nanmean(squeeze(fracFacNeuronCombined(:, nSec, :)))', 'r', 'linew', 3);
    end
    hold off
    xlim([0, 5]);
    ylim([0, 1.5]);
    set(gca, 'Xtick', 1.5:2:5.5, 'XtickLabel', {'A', 'P'});
    ylabel('#factored/#active')
    title(fishListType{2, nExpType})
    fracFacNeuronRatio = squeeze(fracFacNeuronCombined(:, :, 2)./fracFacNeuronCombined(:, :, 1));
    pFracFacNeuronAblaiton(nExpType) = signrank(fracFacNeuronRatio(:, 1), fracFacNeuronRatio(:, 2));
    
    subplot(2, size(fishListType, 2), size(fishListType, 2)+nExpType);
    hold on,
    for nSec = 1:2
        xSeq = repmat([2*nSec-1; 2*nSec], 1, numel(currfishList)*2);
        ySeq = squeeze(factorSpanCombined(:, nSec, :))';
        scatter(xSeq(:), ySeq(:), 'ok');
        plot(xSeq, ySeq, 'k');
        p = signrank(ySeq(1, :), ySeq(2, :));
        if p < pThres
            text(2*nSec-0.5, max(ySeq(:))*1.2, '*');
        end
        scatter([2*nSec-1, 2*nSec], nanmean(squeeze(factorSpanCombined(:, nSec, :)))', 'or', 'linew', 3);
        plot([2*nSec-1, 2*nSec], nanmean(squeeze(factorSpanCombined(:, nSec, :)))', 'r', 'linew', 3);
    end
    hold off
    xlim([0, 5]);
    ylim([-0.1, 1.5])
    set(gca, 'Xtick', 1.5:2:5.5, 'XtickLabel', {'A', 'P'});
    ylabel('normalized AP span')
    title(fishListType{2, nExpType})
    factorSpanRatio = squeeze(factorSpanCombined(:, :, 2)./factorSpanCombined(:, :, 1));
    pFactorSpanAblation(nExpType) = signrank(factorSpanRatio(:, 1), factorSpanRatio(:, 2));
    ratioList = [ratioList; [fracFacNeuronRatio, factorSpanRatio, ones(size(factorSpanRatio, 1), 1)*nExpType] ];
    ratioStats(nExpType, :) = [mean(fracFacNeuronRatio), std(fracFacNeuronRatio)];
    spanStats(nExpType, :) = [mean(factorSpanRatio), std(factorSpanRatio)];
end
catch
    disp('');
end

setPrint(8*3, 6*2, [plotDir '/AblationTypeSummary' tagExt], 'pdf');

mColor = lines(2);
figure, 
errorbar(1:3, ratioStats(:, 1), ratioStats(:, 3), 'color', mColor(1, :), 'linew', 2);
% plot(1:3, ratioStats(:, 1), 'color', mColor(1, :), 'linew', 2);
hold on
errorbar(1:3, ratioStats(:, 2), ratioStats(:, 4), 'color', mColor(2, :), 'linew', 2);
% plot(1:3, ratioStats(:, 2), 'color', mColor(2, :), 'linew', 2);
scatter(ratioList(:, 5)+rand(size(ratioList(:, 5)))*0.2-0.1, ratioList(:, 1), 'filled', 'MarkerFaceColor', mColor(1, :), 'MarkerFaceAlpha', .5);
scatter(ratioList(:, 5)+rand(size(ratioList(:, 5)))*0.2-0.1, ratioList(:, 2), 'filled', 'MarkerFaceColor', mColor(2, :), 'MarkerFaceAlpha', .5);
set(gca ,'XTick', 1:3, 'XTickLabel', {'Cut A', 'Cut M', 'Cut P'});
set(gca, 'YTick', 0:0.4:0.8);
ylabel('Ratio of synchronization level after vs. before cut')
legend('Anterior to cut', 'Posterior to cut');
box off
setPrint(16, 12, [plotDir '/AblationTypeSyncRatio' tagExt], 'pdf');

figure, 
errorbar(1:3, spanStats(:, 1), spanStats(:, 3), 'color', mColor(1, :), 'linew', 2);
% plot(1:3, ratioStats(:, 1), 'color', mColor(1, :), 'linew', 2);
hold on
errorbar(1:3, spanStats(:, 2), spanStats(:, 4), 'color', mColor(2, :), 'linew', 2);
% plot(1:3, ratioStats(:, 2), 'color', mColor(2, :), 'linew', 2);
scatter(ratioList(:, 5)+rand(size(ratioList(:, 5)))*0.2-0.1, ratioList(:, 3), 'filled', 'MarkerFaceColor', mColor(1, :), 'MarkerFaceAlpha', .5);
scatter(ratioList(:, 5)+rand(size(ratioList(:, 5)))*0.2-0.1, ratioList(:, 4), 'filled', 'MarkerFaceColor', mColor(2, :), 'MarkerFaceAlpha', .5);
set(gca ,'XTick', 1:3, 'XTickLabel', {'Cut A', 'Cut M', 'Cut P'});
set(gca, 'YTick', 0:0.4:0.8);
ylabel('Ratio of factor span after vs. before cut')
legend('Anterior to cut', 'Posterior to cut');
box off
setPrint(16, 12, [plotDir '/AblationTypeSpanRatio' tagExt], 'pdf');

close all


