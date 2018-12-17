%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4  (double cut) evaluate change of patterned activity from FA result
%
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%
function Ablation_v_4_2(fishList, tagExt)
addpath('../Func');
setDir;

pThres = 0.05;
nFileList = 25:2:96;


numCells = zeros(numel(fishList), 6, 2);
numActCells = zeros(numel(fishList), 6, 2);
fracFacNeuron = zeros(numel(fishList), 6, 2);
factorSpan = zeros(numel(fishList), 6, 2);



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
        
        nCells = numel(new_x);
        nTimes = size(factorSizes, 1);
        ind = false(nCells, 6);
        ind(:, 1) = new_x>0 & new_x<segAblation_A(1) & new_y<0;
        ind(:, 2) = new_x>0 & new_x<segAblation_A(1) & new_y>0;
        ind(:, 3) = new_x>segAblation_A(2) & new_x<segAblation_P(1) & new_y<0;
        ind(:, 4) = new_x>segAblation_A(2) & new_x<segAblation_P(1) & new_y>0;
        ind(:, 5) = new_x>segAblation_P(2) & new_y<0;
        ind(:, 6) = new_x>segAblation_P(2) & new_y>0;
        
        numActCellts    = nan(nTimes, 6);
        fracFacNeuronts = nan(nTimes, 6);
        factorSpants    = nan(nTimes, 6);
        % method for using all segments
        for nTime = 1:nTimes
            activeTag = activeNeuronMat(:, nTime)>0;
            factorTag = false(numel(activeTag), 1);
            for nFactor = 1:size(factorSizes, 2)
                currentFactor = factorComp(:, nFactor);
                factorTag(currentFactor{nTime}) = 1;
            end
            for nSec = 1:6
                factorTagSec = factorTag & ind(:, nSec);
                activeTagSec = activeTag & ind(:, nSec);
                if sum(factorTagSec)>0 && sum(activeTagSec)>0
                    numActCellts(nTime, nSec) = nanmax(numActCellts(nTime, nSec), sum(activeTagSec));
                    fracFacNeuronts(nTime, nSec) = nanmax(fracFacNeuronts(nTime, nSec), sum(factorTagSec)/sum(activeTagSec));
                    factorSpants(nTime, nSec)   = nanmax(factorSpants(nTime, nSec), range(new_x(factorTagSec))/range(new_x(activeTagSec)));
                end
            end
        end
        %             % alternative method for post ablation: calculate dominate factor
        %             for nFactor = 1:size(factorSizes, 2)
        %                 currentFactor = factorComp(:, nFactor);
        %                 for nTime = 1:nTimes
        %                     activeTag = activeNeuronMat(:, nTime)>0;
        %                     factorTag = false(numel(activeTag), 1);
        %                     factorTag(currentFactor{nTime}) = 1;
        %                     [~, domSec] = max(sum(bsxfun(@and, factorTag, ind)));
        %                     factorTag = factorTag & ind(:, domSec);
        %                     activeTag = activeTag & ind(:, domSec);
        %                     if  sum(factorTag)>0 && numel(domSec)==1
        %                         numActCellts(nTime, domSec) = nanmax(numActCellts(nTime, domSec), sum(activeTag));
        %                         fracFacNeuronts(nTime, domSec) = nanmax(fracFacNeuronts(nTime, domSec), sum(factorTag)/sum(activeTag));
        %                         try
        %                         factorSpants(nTime, domSec)   = nanmax(factorSpants(nTime, domSec), range(new_x(factorTag))/range(new_x(activeTag)));
        %                         catch
        %                             disp('');
        %                         end
        %                     end
        %                 end
        
        numActCells(i, :, nExp) = nanmedian(numActCellts, 1);
        numCells(i, :, nExp) = sum(ind);
        fracFacNeuron(i, :, nExp) = nanmedian(fracFacNeuronts, 1);
        factorSpan(i, :, nExp) = nanmedian(factorSpants, 1);
    end
end
[inv1, inv2] = ind2sub([numel(fishList), 6], find(isnan(fracFacNeuron(:, :, 1)) & isnan(fracFacNeuron(:, :, 2))));
fracFacNeuron(isnan(fracFacNeuron)) = 0;
fracFacNeuron(inv1, inv2, :) = NaN;

[inv1, inv2] = ind2sub([numel(fishList), 6], find(isnan(factorSpan(:, :, 1)) & isnan(factorSpan(:, :, 2))));
factorSpan(isnan(factorSpan)) = 0;
factorSpan(inv1, inv2, :) = NaN;

% figure 1 fracFacNeuron and factorSpan before & after ablation
figure('Position', [0, 0, 400*2, 300]);

fracFacNeuronCombined = [fracFacNeuron(:, [1, 3, 5], :); fracFacNeuron(:, [2, 4, 6], :)];
factorSpanCombined = [factorSpan(:, [1, 3, 5], :); factorSpan(:, [2, 4, 6], :)];
subplot(1, 2, 1);
hold on,
for nSec = 1:3
    xSeq = repmat([2*nSec-1; 2*nSec], 1, numel(fishList)*2);
    ySeq = squeeze(fracFacNeuronCombined(:, nSec, :))';
    scatter(xSeq(:), ySeq(:), 'ok');
    plot(xSeq, ySeq, 'k');
    p = signrank(ySeq(1, :), ySeq(2, :));
    if p<pThres
        text(2*nSec-0.5, max(ySeq(:))*1.2, '*');
        disp(['nSec=' num2str(nSec) ' ,p=' num2str(p)]);
    end
    scatter([2*nSec-1, 2*nSec], nanmean(squeeze(fracFacNeuronCombined(:, nSec, :)))', 'or', 'linew', 3);
    plot([2*nSec-1, 2*nSec], nanmean(squeeze(fracFacNeuronCombined(:, nSec, :)))', 'r', 'linew', 3);
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


fracFacNeuronRatio = squeeze(fracFacNeuronCombined(:, :, 2)./fracFacNeuronCombined(:, :, 1));
fracFacNeuronRatio(fracFacNeuronRatio>1.2) = NaN;
fracFacNeuronRatio(sum(isnan(fracFacNeuronRatio), 2)>0, :) = [];
pFracFacNeuronRatio.AP = signrank(fracFacNeuronRatio(:, 1), fracFacNeuronRatio(:, 3));
pFracFacNeuronRatio.AM = signrank(fracFacNeuronRatio(:, 1), fracFacNeuronRatio(:, 2));
pFracFacNeuronRatio.MP = signrank(fracFacNeuronRatio(:, 2), fracFacNeuronRatio(:, 3));
pFracFacNeuronRatio.M_AP = signrank(fracFacNeuronRatio(:, 2), nanmax(fracFacNeuronRatio(:, [1, 3]), [], 2));

factorSpanRatio = squeeze(factorSpanCombined(:, :, 2)./factorSpanCombined(:, :, 1));
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
title({['pAM = ', num2str(pFracFacNeuronRatio.AM)], ['pAP = ', num2str(pFracFacNeuronRatio.AP)], ['pMP = ', num2str(pFracFacNeuronRatio.MP)]});
set(gca, 'XTick', 1:3, 'XTickLabel', {'A', 'M', 'P'});
ylabel({'ratio of #factored/#active', 'after vs. before'});
subplot(1, 2, 2)
hold on
plot(1:3, factorSpanRatio, 'ko-');
plot(1:3, nanmean(factorSpanRatio), 'ro-', 'linew', 2);
title({['pAM = ', num2str(pFactorSpanAblation.AM)], ['pAP = ', num2str(pFactorSpanAblation.AP)], ['pMP = ', num2str(pFactorSpanAblation.MP)]});
set(gca, 'XTick', 1:3, 'XTickLabel', {'A', 'M', 'P'});
ylabel({'ratio of factor span', 'after vs. before'});
setPrint(16, 6, [plotDir '/AblationTypeSyncRatio' tagExt], 'pdf');



