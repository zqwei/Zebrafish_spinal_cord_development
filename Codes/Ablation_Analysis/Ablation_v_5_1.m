%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5  (double cut) matrix plot for individual fish
%
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%

function Ablation_v_5_1(fishList)
addpath('../Func');
setDir;


fsSpanThresFull = 0.5; % factor span exceeds this size, call global
fsSpanThresPart = 0;% factor span under this size, call non-existing

factorSpanAll = nan(numel(fishList)*2, 3);
lifeTimeGlobal = nan(numel(fishList)*2, 3);

for i = 1:numel(fishList);
    nFile = 24+fishList(i)*2;
    fileDirName  = fileDirNames{nFile};
    fileName          = fileNames{nFile}; 

    load([fileDirName '/', 'profile.mat'], 'segAblation_A', 'segAblation_P');
    load([tempDatDir, fileName, '.mat'], 'new_x', 'new_y', 'activeNeuronMat'); 
    load([tempDatDir, 'LONOLoading_' fileName, '_v2.mat'], 'factorComp', 'factorSizes');
    

    invalidIndex = (new_x>segAblation_A(1) & new_x<=segAblation_A(2)) | (new_x>segAblation_P(1) & new_x<=segAblation_P(2));
    
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
    factorSpanAll(2*i-1, :) = nanmedian(factorSpan(:, 1:3));
    factorSpanAll(2*i, :) = nanmedian(factorSpan(:, 4:6));
%     fracActNeuronAll(2*i-1, :) = nanmedian(fracActNeuron(:, 1:3));
%     fracActNeuronAll(2*i, :) = nanmedian(fracActNeuron(:, 4:6));
    lifeTimeGlobal(2*i-1, :) = sum(factorSpan(:, 1:3)>fsSpanThresFull)/size(factorSpan, 1);
    lifeTimeGlobal(2*i, :) = sum(factorSpan(:, 4:6)>fsSpanThresFull)/size(factorSpan, 1);
end

factorSpanAll(isnan(factorSpanAll)) = 0;


m = factorSpanAll;
m(m>fsSpanThresFull) = 2;
m(m<fsSpanThresFull & m>fsSpanThresPart) = 1;
m(m<fsSpanThresPart) = 0;
im = plotMatrixType(m', 1);
imwrite(im, [plotDir, 'DoubleAblationMatrix_HalfFish.tif']);

m = max(factorSpanAll(1:2:end, :), factorSpanAll(2:2:end, :));
% m = (factorSpanAll(1:2:end, :)+ factorSpanAll(2:2:end, :))/2;
m(m>fsSpanThresFull) = 2;
m(m<fsSpanThresFull & m>fsSpanThresPart) = 1;
m(m<fsSpanThresPart) = 0;
im = plotMatrixType(m', 0);
imwrite(im, [plotDir, 'DoubleAblationMatrix_WholeFish.tif']);


lifeTimeGlobal = max(lifeTimeGlobal(1:2:end, :), lifeTimeGlobal(2:2:end, :));
% lifeTimeGlobal = lifeTimeGlobal(1:2:end, :)+ lifeTimeGlobal(2:2:end, :);
lifeTime_M = lifeTimeGlobal(:, 2);
lifeTime_AP = max(lifeTimeGlobal(:, [1, 3]), [], 2);
pM_AP = signrank(lifeTime_AP, lifeTime_M);

% figure, imagesc(lifeTimeGlobal);
% set(gca, 'XTick', 1:3, 'XTickLabel', {'A', 'M', 'P'});
% ylabel('fish index');
% h = colorbar;
% ylabel(h, 'global factor life time (min)');
% title(['AP vs. M, p=' num2str(pM_AP, '%.4f')]);
% setPrint(8, 6, [plotDir '/DoubleCutFactorLife'], 'pdf');

figure, plot([zeros(numel(lifeTime_M), 1), ones(numel(lifeTime_M), 1)]', [lifeTime_M, lifeTime_AP]', 'ok-');
xlim([-0.5 1.5]);
ylim([0, 1]);
set(gca, 'XTick', [0 1], 'XTickLabel', {'M Factor', 'A/P Factor'});
ylabel('Global factor life time');
title(['AP vs. M, p=' num2str(pM_AP, '%.4f')]);
setPrint(8, 6, [plotDir '/DoubleCutFactorLife'], 'pdf');