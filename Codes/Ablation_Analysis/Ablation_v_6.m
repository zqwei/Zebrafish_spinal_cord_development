%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6  (double cut) region size vs. factor size
% matrix plot with region size info
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%

function Ablation_v_6(fishList, tagExt)
addpath('../Func');
setDir;

fsSpanThresFull = 0.5; % factor span exceeds this size, call global
fsSpanThresPart = 0.1;% factor span under this size, call non-existing

try
factorSpanAll  = nan(numel(fishList)*2, 3);
fracActNeuronAll = nan(numel(fishList)*2, 3);
regionSpanAbs     = nan(numel(fishList)*2, 3);
regionSpan     = nan(numel(fishList)*2, 3);

for i = 1:numel(fishList)
    nFile = 24+fishList(i)*2;
    fileDirName  = fileDirNames{nFile};
    fileName          = fileNames{nFile}; 

    load([fileDirName '/', 'profile.mat'], 'segAblation_A', 'segAblation_P', 'segOffset');
    load([tempDatDir, fileName, '.mat'], 'new_x', 'new_y', 'activeNeuronMat'); 
    load([tempDatDir, 'LONOLoading_' fileName, '_v2.mat'], 'factorComp', 'factorSizes');
    

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
    factorSpanAll(2*i-1, :) = nanmedian(factorSpan(:, 1:3));
    factorSpanAll(2*i, :) = nanmedian(factorSpan(:, 4:6));
    fracActNeuronAll(2*i-1, :) = nanmedian(fracActNeuron(:, 1:3));
    fracActNeuronAll(2*i, :) = nanmedian(fracActNeuron(:, 4:6));
    
% region span is defined by AP span of neurons before ablation
    nFile = 23+fishList(i)*2;
    fileName          = fileNames{nFile}; 

    load([tempDatDir, fileName, '.mat'], 'new_x', 'new_y', 'activeNeuronMat'); 
    load([tempDatDir, 'LONOLoading_' fileName, '_v2.mat'], 'factorComp', 'factorSizes');
    pEnd = max(new_x(cell2mat(factorComp(:))));
    aEnd = min(new_x(cell2mat(factorComp(:))));
    regionSpan(2*i-1:2*i, 1) = segAblation_A(1) - aEnd;
    regionSpan(2*i-1:2*i, 2) = segAblation_P(1) - segAblation_A(2);
    regionSpan(2*i-1:2*i, 3) = pEnd - segAblation_P(2);
   
    regionSpanAbs(2*i-1:2*i, 1) = segAblation_A(1) - 0 + segOffset;
    regionSpanAbs(2*i-1:2*i, 2) = segAblation_P(1) - segAblation_A(2);
    regionSpanAbs(2*i-1:2*i, 3) = 14 - segAblation_P(2) - segOffset;
end

catch
    disp('');
end
factorSpanAll(isnan(factorSpanAll)) = 0;
fracActNeuronAll(isnan(fracActNeuronAll)) = 0;
ablationType = repmat(1:3, numel(fishList)*2, 1);
colorMap = lines(3);

% figure,
% subplot(1, 2, 1)
% scatter(regionSpanAbs(:)+ (rand(numel(regionSpanAbs), 1)-0.5)*0.6, fracActNeuronAll(:),20, colorMap(ablationType(:),:), 'filled');
% hold on
% for xs = 0.5:1:9
%     plot([xs, xs], [0, 0.8], 'k--');
% end
% hold off
% xlabel('Region span (seg)')
% ylabel('#factored/#active');
% 
% 
% subplot(1, 2, 2);
% scatter(regionSpanAbs(:)+ (rand(numel(regionSpanAbs), 1)-0.5) *0.6, factorSpanAll(:), 20, colorMap(ablationType(:),:), 'filled');
% hold on
% for xs = 0.5:1:9
%     plot([xs, xs], [0, 1], 'k--');
% end
% hold off
% xlabel('Region span (seg)')
% ylabel('Relative factor span');
% setPrint(32, 12, [plotDir '/NetworkSize' tagExt], 'pdf');

% m = factorSpanAll;
% m(m>fsSpanThresFull) = 2;
% m(m<fsSpanThresFull & m>fsSpanThresPart) = 1;
% m(m<fsSpanThresPart) = 0;
% im = plotMatrixLocation(m', regionSpanAbs', 1);
% imwrite(im, [plotDir, 'DoubleAblationMatrix_HalfFish' tagExt '.tif']);
% 
% m = max(factorSpanAll(1:2:end, :), factorSpanAll(2:2:end, :));
% regionSpanAbs = regionSpanAbs(1:2:end, :);
% m(m>fsSpanThresFull) = 2;
% m(m<fsSpanThresFull & m>fsSpanThresPart) = 1;
% m(m<fsSpanThresPart) = 0;
% im = plotMatrixLocation(m', regionSpanAbs', 0);
% imwrite(im, [plotDir, 'DoubleAblationMatrix_WholeFish' tagExt '.tif']);
% 
% figure, imagesc(regionSpan(1:2:end, :)'>1);

regionSpanAbs = regionSpanAbs(:);
factorSpanAll = factorSpanAll(:);
ablationType  = ablationType(:);
regionSpanAbsSeq = min(regionSpanAbs):max(regionSpanAbs);
% bins = 0:0.2:1.2;
bins = 0:0.05:1;
densityMap = nan(numel(regionSpanAbsSeq), numel(bins), 3);

figure,
for i = 1:numel(regionSpanAbsSeq)
    subplot(numel(regionSpanAbsSeq), 1, i)
    hold on
    for type = 1:3
        if sum(regionSpanAbs==regionSpanAbsSeq(i) & ablationType==type)>0
%             %option 1: no smoothing raw data
%             f = histcounts(factorSpanAll(regionSpanAbs==regionSpanAbsSeq(i) & ablationType==type), bins);
%             option 2: smoothing by ksdensity
            [f, ~] = ksdensity(factorSpanAll(regionSpanAbs==regionSpanAbsSeq(i) & ablationType==type), bins, 'bandwidth', 0.1);
            plot(bins(1:end), f, 'Color', colorMap(type, :));
            densityMap(i, :, type) = f;
        end
    end
%     xlabel('Normalized AP span (seg)');
%     ylabel({'density'});
%     title(['Region size = ' num2str(regionSpanAbsSeq(i)) ' segs']);
    set(gca, 'XTick', 0:0.2:1);
    hold off
end
setPrint(16, 24, [plotDir '/NetworkSizeByRegionSize' tagExt], 'pdf');


typeText = {'A', 'M', 'P'};
figure, 
for type = 1:3
    subplot(3, 1, type)
    imagesc(regionSpanAbsSeq, bins, densityMap(:, :, type)');
    h = colorbar;
    title(['region type ' typeText{type}])
    xlabel('Region size (seg)')
    ylabel({'Normalized', 'AP span (seg)'})
    set(gca, 'YDir', 'Normal')
    ylabel(h, 'count')
end
setPrint(16, 12, [plotDir '/NetworkSizeByRegionType' tagExt], 'pdf');

