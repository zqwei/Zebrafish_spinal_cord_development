%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4  (double cut) evaluate change of patterned activity from FA result
% use consitent definition as Ablation_v_5_1
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%
function Ablation_v_4_3(fishList)
addpath('../Func');
setDir;

pThres = 0.05;
nFileList = 25:2:96;


factorSpanAll  = nan(numel(fishList)*2, 3);
fracActNeuronAll = nan(numel(fishList)*2, 3);

tagExt = '_Double_FA_JointWin_DefinedSeg';


for i = 1:numel(fishList)
    nFish = fishList(i);
    for nExp = 1:2
        nFile = nFileList(nFish) + nExp - 1;
        fileDirName   = fileDirNames{nFile}; %#ok<*USENS>
        fileName      = fileNames{nFile};
        
        dirImageData  = [fileDirName '/'];
        load([dirImageData, 'profile.mat'], 'segAblation_A', 'segAblation_P');
        load([tempDatDir, 'LONOLoading_' fileName, '_v2.mat'], 'factorComp', 'factorSizes');
        load([tempDatDir, fileName, '.mat'], 'dff', 'activeNeuronMat', 'new_x', 'new_y');
        
        invalidIndex = (new_x<=0 & new_x>segAblation_A(1) & new_x<=segAblation_A(2)) | (new_x>segAblation_P(1) & new_x<=segAblation_P(2));
        
        factorSpan    = nan(size(factorSizes, 1), 6);
        fracActNeuron = nan(size(factorSizes, 1), 6);
        for nFactor = 1:size(factorSizes, 2)
            currentFactor = factorComp(:, nFactor);
            for nTime = 1:numel(currentFactor)
                validIndex = currentFactor{nTime};
                validIndex(invalidIndex(validIndex)) = [];
                xFac = new_x(validIndex);
                if isempty(xFac)
                    continue;
                end
                sideFac = unique(new_y(currentFactor{nTime})>0);
                [voteAMP, segEdges] = histcounts(xFac, [-Inf, mean(segAblation_A), mean(segAblation_P), Inf]);
                [~, domSec] = max(voteAMP);
                xFac     = xFac(xFac>segEdges(domSec) & xFac<=segEdges(domSec+1));
                currFrac = numel(xFac)/sum(activeNeuronMat(new_x>segEdges(domSec) & new_x<=segEdges(domSec+1) & ~invalidIndex, nTime));
                currFS = range(xFac);
                % normalize by section length
                xLen = range(new_x(new_x>segEdges(domSec) & new_x<=segEdges(domSec+1) & ~invalidIndex & activeNeuronMat(:, nTime)));
                frac = sum(sum(new_x>segEdges(domSec) & new_x<=segEdges(domSec+1) & ~invalidIndex & activeNeuronMat(:, nTime)));
                currFS = currFS/xLen;
                factorSpan(nTime, domSec + sideFac*3) = nanmax(factorSpan(nTime, domSec + sideFac*3), currFS);
                fracActNeuron(nTime, domSec + sideFac*3) = nanmax(fracActNeuron(nTime, domSec + sideFac*3), currFrac);
            end
        end
        %     figure, imagesc(factorSpan, [0, 2]); colorbar
        factorSpanAll(2*i-1, :, nExp) = nanmedian(factorSpan(:, 1:3));
        factorSpanAll(2*i, :, nExp) = nanmedian(factorSpan(:, 4:6));
        fracActNeuronAll(2*i-1, :, nExp) = nanmedian(fracActNeuron(:, 1:3));
        fracActNeuronAll(2*i, :, nExp) = nanmedian(fracActNeuron(:, 4:6));
    end
end

fracActNeuronAll(isnan(fracActNeuronAll)) = 0;
factorSpanAll(isnan(factorSpanAll)) = 0;

% figure 1 fracFacNeuron and factorSpan before & after ablation
figure('Position', [0, 0, 400*2, 300]);

subplot(1, 2, 1);
hold on,
for nSec = 1:3
    xSeq = repmat([2*nSec-1; 2*nSec], 1, numel(fishList)*2);
    ySeq = squeeze(fracActNeuronAll(:, nSec, :))';
    scatter(xSeq(:), ySeq(:), 'ok');
    plot(xSeq, ySeq, 'k');
    p = signrank(ySeq(1, :), ySeq(2, :));
    if p<pThres
        text(2*nSec-0.5, max(ySeq(:))*1.2, '*');
    end
    scatter([2*nSec-1, 2*nSec], nanmean(squeeze(fracActNeuronAll(:, nSec, :)))', 'or', 'linew', 3);
    plot([2*nSec-1, 2*nSec], nanmean(squeeze(fracActNeuronAll(:, nSec, :)))', 'r', 'linew', 3);
end
hold off
xlim([0, 7]);
ylim([-0.2, 1.5]);
set(gca, 'Xtick', 1.5:2:5.5, 'XtickLabel', {'A', 'M', 'P'});
ylabel('#factored/#active')


subplot(1, 2, 2);
hold on,
for nSec = 1:3
    xSeq = repmat([2*nSec-1; 2*nSec], 1, numel(fishList)*2);
    ySeq = squeeze(factorSpanAll(:, nSec, :))';
    scatter(xSeq(:), ySeq(:), 'ok');
    plot(xSeq, ySeq, 'k');
    p = signrank(ySeq(1, :), ySeq(2, :));
    if p < pThres
        text(2*nSec-0.5, max(ySeq(:))*1.2, '*');
    end
    scatter([2*nSec-1, 2*nSec], nanmean(squeeze(factorSpanAll(:, nSec, :)))', 'or', 'linew', 3);
    plot([2*nSec-1, 2*nSec], nanmean(squeeze(factorSpanAll(:, nSec, :)))', 'r', 'linew', 3);
end
hold off
xlim([0, 7]);
ylim([-0.1, 1.5])
set(gca, 'Xtick', 1.5:2:5.5, 'XtickLabel', {'A', 'M', 'P'});
ylabel('normalized AP span')

setPrint(25, 6, [plotDir '/AblationTypeSummary' tagExt], 'pdf');

% figure showing synchoronization level change
% figure('Position', [0, 0, 400*2, 300]);
% subplot(1, 2, 1);
% imagesc((factorSpanCombined(1:8, :, 1)+factorSpanCombined(9:end, :, 1))/2, [0, 1]);
% set(gca, 'Xtick', [1 2 3], 'XtickLabel', {'A', 'M', 'P'});
% colorbar
% title('before ablation')
% ylabel('fish number');
% subplot(1, 2, 2);
% imagesc((factorSpanCombined(1:8, :, 2)+factorSpanCombined(9:end, :, 2))/2, [0, 1]);
% set(gca, 'Xtick', [1 2 3], 'XtickLabel', {'A', 'M', 'P'});
% ylabel('fish number');
% colorbar
% title('after ablation');


fracFacNeuronRatio = squeeze(fracActNeuronAll(:, :, 2)./fracActNeuronAll(:, :, 1));
fracFacNeuronRatio(fracFacNeuronRatio>1.2) = NaN;
fracFacNeuronRatio(sum(isnan(fracFacNeuronRatio), 2)>0, :) = [];
pFracFacNeuronRatio.AP = signrank(fracFacNeuronRatio(:, 1), fracFacNeuronRatio(:, 3));
pFracFacNeuronRatio.AM = signrank(fracFacNeuronRatio(:, 1), fracFacNeuronRatio(:, 2));
pFracFacNeuronRatio.MP = signrank(fracFacNeuronRatio(:, 2), fracFacNeuronRatio(:, 3));
pFracFacNeuronRatio.M_AP = signrank(fracFacNeuronRatio(:, 2), nanmax(fracFacNeuronRatio(:, [1, 3]), [], 2));

factorSpanRatio = squeeze(factorSpanAll(:, :, 2)./factorSpanAll(:, :, 1));
factorSpanRatio(factorSpanRatio>1.2) = NaN;
factorSpanRatio(sum(isnan(factorSpanRatio), 2)>0, :) = [];
pFactorSpanAblation.AP = signrank(factorSpanRatio(:, 1), factorSpanRatio(:, 3));
pFactorSpanAblation.AM = signrank(factorSpanRatio(:, 1), factorSpanRatio(:, 2));
pFactorSpanAblation.MP = signrank(factorSpanRatio(:, 2), factorSpanRatio(:, 3));
pFactorSpanAblation.M_AP = signrank(factorSpanRatio(:, 2), nanmax(factorSpanRatio(:, [1, 3]), [], 2));

figure,
subplot(1, 2, 1)
hold on
plot(1:3, fracFacNeuronRatio, 'ko-');
plot(1:3, nanmean(fracFacNeuronRatio), 'ro-', 'linew', 2);
title(['p-value M vs. AP = ', num2str(pFracFacNeuronRatio.M_AP)]);
set(gca, 'XTick', 1:3, 'XTickLabel', {'A', 'M', 'P'});
ylabel('change of factored neuron ratio after vs. before');
subplot(1, 2, 2)
hold on
plot(1:3, factorSpanRatio, 'ko-');
plot(1:3, nanmean(factorSpanRatio), 'ro-', 'linew', 2);
title(['p-value M vs. AP = ', num2str(pFactorSpanAblation.M_AP)]);
set(gca, 'XTick', 1:3, 'XTickLabel', {'A', 'M', 'P'});
ylabel('change of factor span after vs. before');
setPrint(16, 6, [plotDir '/AblationTypeSyncRatio' tagExt], 'pdf');



