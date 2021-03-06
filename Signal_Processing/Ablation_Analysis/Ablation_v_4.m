%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4  (double cut) evaluate change of patterned activity from FA result
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

numCells = zeros(numel(fishList), 4, 2);
numActCells = zeros(numel(fishList), 4, 2);
fracFacNeuron = zeros(numel(fishList), 4, 2);
factorSpan = zeros(numel(fishList), 4, 2);

tagExt = '_Double_FA';


for i = 1:numel(fishList)
    nFish = fishList(i);
    for nExp = 1:2
        nFile = nFileList(nFish) + nExp - 1;
        fileDirName   = fileDirNames{nFile}; %#ok<*USENS>
        fileName      = fileNames{nFile};
        
        dirImageData  = [fileDirName '/'];
        load([dirImageData, 'profile.mat'], 'segAblation_A', 'segAblation_P');
        load([tempDatDir, 'FALONO_Average_' fileName, '.mat'], 'LMat');
        load([tempDatDir, fileName, '.mat'], 'dff', 'activeNeuronMat', 'new_x', 'new_y');
        
        ind{1} = new_x<segAblation_A(1) & new_y<0;
        ind{2} = new_x<segAblation_A(1) & new_y>0;
        ind{3} = new_x>segAblation_A(2) & new_x<segAblation_P(1) & new_y<0;
        ind{4} = new_x>segAblation_A(2) & new_x<segAblation_P(1) & new_y>0;
        ind{5} = new_x>segAblation_P(2) & new_y<0;
        ind{6} = new_x>segAblation_P(2) & new_y>0;
        activeTag = sum(activeNeuronMat, 2)>0;
        factorTag = sum(LMat, 2) >0;
        %         factorTag = sum(double(LMat>0) .* repmat(1:size(LMat, 2), size(LMat, 1), 1), 2);
        
        % 1. fraction of factored neruons in the 4 sectors
        for nSec = 1:6
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
    end
    scatter([2*nSec-1, 2*nSec], nanmean(squeeze(fracFacNeuronCombined(:, nSec, :)))', 'or', 'linew', 3);
    plot([2*nSec-1, 2*nSec], nanmean(squeeze(fracFacNeuronCombined(:, nSec, :)))', 'r', 'linew', 3);
end
hold off
xlim([0, 7]);
ylim([-0.2, 1.5]);
set(gca, 'Xtick', 1.5:2:5.5, 'XtickLabel', {'A', 'M', 'P'});
ylabel('#factored/#active')
fracFacNeuronRatio = squeeze(fracFacNeuronCombined(:, :, 1) - fracFacNeuronCombined(:, :, 2));

pFracFacNeuronRatio.AP = signrank(fracFacNeuronRatio(:, 1), fracFacNeuronRatio(:, 3));
pFracFacNeuronRatio.AM = signrank(fracFacNeuronRatio(:, 1), fracFacNeuronRatio(:, 2));
pFracFacNeuronRatio.MP = signrank(fracFacNeuronRatio(:, 2), fracFacNeuronRatio(:, 3));

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
factorSpanRatio = squeeze(factorSpanCombined(:, :, 1) - factorSpanCombined(:, :, 2));
pFactorSpanAblation.AP = signrank(factorSpanRatio(:, 1), factorSpanRatio(:, 3));
pFactorSpanAblation.AM = signrank(factorSpanRatio(:, 1), factorSpanRatio(:, 2));
pFactorSpanAblation.MP = signrank(factorSpanRatio(:, 2), factorSpanRatio(:, 3));

setPrint(8*3, 6*2, [plotDir '/AblationTypeSummary' tagExt], 'pdf');

figure('Position', [0, 0, 400*2, 300]);
subplot(1, 2, 1);
imagesc((factorSpanCombined(1:9, :, 1)+factorSpanCombined(10:end, :, 1))/2, [0, 1]);
set(gca, 'Xtick', [1 2 3], 'XtickLabel', {'A', 'M', 'P'});
colorbar
title('before ablation')
ylabel('fish number');
subplot(1, 2, 2);
imagesc((factorSpanCombined(1:9, :, 2)+factorSpanCombined(10:end, :, 2))/2, [0, 1]);
set(gca, 'Xtick', [1 2 3], 'XtickLabel', {'A', 'M', 'P'});
ylabel('fish number');
colorbar
title('after ablation');

% close all

% cross table for independence test
factorTag = fracFacNeuronCombined(:, :, 2)>0;
[hFisher.AM, pFisher.AM , ~] = fishertest(crosstab(factorTag(:, 1), factorTag(:, 2)));
[hFisher.AP, pFisher.AP , ~] = fishertest(crosstab(factorTag(:, 1), factorTag(:, 3)));
[hFisher.PM, pFisher.PM , ~] = fishertest(crosstab(factorTag(:, 2), factorTag(:, 3)));

[hFisher.EitherEdge, pFisher.EitherEdge , ~] = fishertest(crosstab(factorTag(:, 2), factorTag(:, 3)|factorTag(:, 1)));
[hFisher.BothEdge, pFisher.BothEdge , ~] = fishertest(crosstab(factorTag(:, 2), factorTag(:, 3)&factorTag(:, 1)));

n_boot = 10000;
bootMean = bootstrp(n_boot, @mean, factorTag);

p_AM = mean(bootMean(:, 1) < bootMean(:, 2));
p_PM = mean(bootMean(:, 3) < bootMean(:, 2));
p_AP = mean(bootMean(:, 1) < bootMean(:, 3));

