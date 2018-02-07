%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.  evaluate change of patterned activity from FA result
%
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%

addpath('../Func');
setDir;

pThres = 0.05;
nFileList = 25:2:58;

fishListCutA = [4, 2, 3, 7, 1]; % anterior cut
fishListCutM = [12, 9, 10, 8]; % middle cut
fishListCutP = [16, 17, 13, 15]; % posterior cut
% sectorNames = {'AL', 'AR', 'PL', 'PR'};
sectorNames = {'Anterior', 'Posterior'};
tagExt = [];

fishList = [fishListCutA, fishListCutM, fishListCutP]; % anterior cut

fracPatternNeuron = zeros(numel(fishList), 2, 2);
factorSizeMean = zeros(numel(fishList), 2, 2);
factorSizeStd = zeros(numel(fishList), 2, 2);
factorSpanMean = zeros(numel(fishList), 2, 2);
factorSpanStd = zeros(numel(fishList), 2, 2);


for i = 1:numel(fishList)
    nFish = fishList(i);
    for nExp = 1:2
        nFile = nFileList(nFish) + nExp - 1;
        fileDirName   = fileDirNames{nFile}; %#ok<*USENS>
        fileName      = fileNames{nFile};
        
        dirImageData  = [fileDirName '/'];
        load([dirImageData, 'profile.mat'], 'segAblation');
        load([tempDatDir, fileName, '.mat'], 'dff', 'activeNeuronMat', 'new_x', 'new_y', 'new_z', 'slicedIndex');
        load([tempDatDir, 'FALONO_Average_', fileName, '.mat'], 'LMat', 'PsiMat');
        
        if nExp == 1
            activeTag = activeNeuronMat > 0;
        else
            activeTag = sum(activeNeuronMat, 2)>0;
        end
        factorTag = sum(LMat, 2)>0;
        
%         ind{1} = new_x<segAblation(1) & new_y<0;
%         ind{2} = new_x<segAblation(1) & new_y>0;
%         ind{3} = new_x>segAblation(2) & new_y<0;
%         ind{4} = new_x>segAblation(2) & new_y>0;
        ind{1} = new_x<segAblation(1);
        ind{2} = new_x>segAblation(2);
        
        
        for nSec = 1:2
            % 1. number of factored neurons/active neurons
            fracPatternNeuron(i, nSec, nExp) = sum(factorTag & ind{nSec})/sum(activeTag & ind{nSec});
            % 2. Average and standard devation of 1)factor size, 2)factor span
            LMatSec = LMat(ind{nSec}, :)>0;
            xSeq = new_x(ind{nSec}, :);
            LMatSec(sum(LMatSec, 2)==0, :) = [];
            factorSizeMean(i, nSec, nExp) = mean(sum(LMatSec, 1)/sum(activeTag & ind{nSec}));
            factorSizeStd(i, nSec, nExp) = std(sum(LMatSec, 1)/sum(activeTag & ind{nSec}));
            factorSpans = [];
            for nFactor = 1:size(LMatSec, 2)
                xSeqFactor = xSeq(LMatSec(:, nFactor));
                factorSpans = [factorSpans, (max(xSeqFactor)-min(xSeqFactor))/(max(xSeq)-min(xSeq))];
            end
            if isempty(factorSpans)
                factorSpans = 0;
            end
            factorSpanMean(i, nSec, nExp) = mean(factorSpans);
            factorSpanStd(i, nSec, nExp) = std(factorSpans);
        end
    end
end

fishListType = cell(numel(fishList), 1);
fishListType(1:numel(fishListCutA)) = {'Cut A'};
fishListType(numel(fishListCutA)+1:numel(fishListCutA)+numel(fishListCutM)) = {'Cut M'};
fishListType(numel(fishListCutA)+numel(fishListCutM)+1:numel(fishList)) = {'Cut P'};

figure('Position', [0, 0, 400*2, 300]);
for nSec = 1:2
    subplot(1, 2, nSec);
    hold on
    errorbar(1:numel(fishList), squeeze(factorSpanMean(:, nSec, 1)), squeeze(factorSpanStd(:, nSec, 1)), 'b');
    errorbar(1:numel(fishList), squeeze(factorSpanMean(:, nSec, 2)), squeeze(factorSpanStd(:, nSec, 2)), 'r');
    ylabel('span of the factors (segments)');
    xlabel('fish number (ablation site anterior -> posterior)')
    title([sectorNames{nSec}]);
    xlim([0 numel(fishList)+1]);
    ylim([-0.2, 1.5]);
    legend({'before ablation', 'after ablation'});
end
export_fig([plotDir, 'AblationFA_FactorSpan.pdf']);

figure('Position', [0, 0, 400*2, 300]);
for nSec = 1:2
    subplot(1, 2, nSec);
    hold on
    errorbar(1:numel(fishList), squeeze(factorSizeMean(:, nSec, 1)), squeeze(factorSizeStd(:, nSec, 1)), 'b');
    errorbar(1:numel(fishList), squeeze(factorSizeMean(:, nSec, 2)), squeeze(factorSizeStd(:, nSec, 2)), 'r');
    title([sectorNames{nSec}]);
    legend({'before ablation', 'after ablation'});
    xlim([0 numel(fishList)+1]);
    ylim([-0.1, 0.8]);
    ylabel('fraction of factored neurons');
    xlabel('fish number (ablation site anterior -> posterior)');
end
export_fig([plotDir, 'AblationFA_FracFactoredNeurons.pdf']);
close all;