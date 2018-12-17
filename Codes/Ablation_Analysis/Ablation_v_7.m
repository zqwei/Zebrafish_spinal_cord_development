%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6  (single & double cut) factor size (sync level) wrt region size and location
% matrix plot with region size and location
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%

function Ablation_v_7(fishListDouble, fishListSingle)
addpath('../Func');
setDir;



factorSpanAll    = [];
fracActNeuronAll = [];
regionStartAbs   = [];
regionSpanAll    = [];
fishIDAll        = [];
regionTypeAll    = [];%AL,AR,ML,MR,PL,PR: 1-6

lastSeg = 14;
for i = 1:numel(fishListDouble)
    nFile = 24+fishListDouble(i)*2;
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
    factorSpanAll    = [factorSpanAll; nanmedian(factorSpan(:, 1:3))'];
    factorSpanAll    = [factorSpanAll; nanmedian(factorSpan(:, 4:6))'];
    fracActNeuronAll = [fracActNeuronAll; nanmedian(fracActNeuron(:, 1:3))'];
    fracActNeuronAll = [fracActNeuronAll; nanmedian(fracActNeuron(:, 4:6))'];
    
    regionSpan       = [segAblation_A(1) - 0 + segOffset; segAblation_P(1) - segAblation_A(2); lastSeg - segAblation_P(2) - segOffset];
    regionSpanAll    = [regionSpanAll; regionSpan; regionSpan];
    
    regionStart      = [1; segAblation_A(2) + segOffset; segAblation_P(1) + segOffset];
%     regionStart      = round([mean([0, segAblation_A(2) + segOffset]); mean([segAblation_A(2) + segOffset, segAblation_P(1) + segOffset]); mean([segAblation_P(1) + segOffset, lastSeg])]);
    regionStartAbs   = [regionStartAbs; regionStart; regionStart];
    
    fishIDAll        = [fishIDAll; repmat(fishListDouble(i), 6, 1)];
    regionTypeAll    = [regionTypeAll; (1:6)'];
end


for i = 1:numel(fishListSingle)
    nFile = 24+fishListSingle(i)*2;
    fileDirName   = fileDirNames{nFile}; %#ok<*USENS>
    fileName      = fileNames{nFile};
    
    dirImageData  = [fileDirName '/'];
    load([dirImageData, 'profile.mat'], 'segAblation', 'segOffset');
    load([tempDatDir, fileName, '.mat'], 'dff', 'activeNeuronMat', 'new_x', 'new_y', 'new_z', 'slicedIndex');
    load([tempDatDir, 'FALONO_Average_', fileName, '.mat'], 'LMat');
    
    activeTag = sum(activeNeuronMat, 2)>0;
    factorTag = sum(LMat, 2)>0;
    
    ind{1} = new_x<segAblation(1) & new_y<0;
    ind{2} = new_x<segAblation(1) & new_y>0;
    ind{3} = new_x>segAblation(2) & new_y<0;
    ind{4} = new_x>segAblation(2) & new_y>0;
    
    
    for nSec = 1:4
        % 1. number of factored neurons/active neurons
        fracPatternNeuron(i, nSec) = sum(factorTag & ind{nSec})/sum(activeTag & ind{nSec});
        % 2. 1)factor size, 2)factor span
        LMatSec = LMat(ind{nSec}, :)>0;
        xSeq = new_x(ind{nSec}, :);
        LMatSec(:, sum(LMatSec, 1)==0) = [];
        factorSizes = sum(LMatSec, 1)/sum(activeTag & ind{nSec});
        factorSpans = [];
        for nFactor = 1:size(LMatSec, 2)
            xSeqFactor = xSeq(LMatSec(:, nFactor));
            factorSpans = [factorSpans, (max(xSeqFactor)-min(xSeqFactor))/(max(xSeq)-min(xSeq))];
            %             factorSpans = [factorSpans, max(xSeqFactor)-min(xSeqFactor)];
        end
        if isempty(factorSpans)
            factorSpans = 0;
            factorSizes = 0;
        end
        factorSpanAll = [factorSpanAll;max(factorSpans)];
    end
    regionSpan     = [segAblation(1)+segOffset; segAblation(1)+segOffset;lastSeg - segAblation(2) - segOffset;lastSeg - segAblation(2) - segOffset];
    regionSpanAll  = [regionSpanAll; regionSpan];
    
    regionStart    = [1; 1; segAblation(2) + segOffset;segAblation(2) + segOffset];
%     regionStart    = [mean([0, segAblation(2) + segOffset]); mean([0, segAblation(2) + segOffset]); mean([segAblation(2) + segOffset, lastSeg]);mean([segAblation(2) + segOffset, lastSeg]);];
    regionStartAbs = [regionStartAbs; regionStart];
    
    fishIDAll        = [fishIDAll; repmat(fishListSingle(i), 4, 1)];
    regionTypeAll    = [regionTypeAll; 1;1;5;6];
end

factorSpanAll(isnan(factorSpanAll)) = 0;
fracActNeuronAll(isnan(fracActNeuronAll)) = 0;

binsStart      = min(regionStartAbs):max(regionStartAbs);
binsRegionSpan = min(regionSpanAll):max(regionSpanAll);
fsBins = nan(numel(binsStart), numel(binsRegionSpan));
for i = 1:numel(binsStart)
    for j = 1:numel(binsRegionSpan)
        if sum(regionStartAbs==binsStart(i) & regionSpanAll==binsRegionSpan(j))>0
            fsBins(i, j) = nanmedian(factorSpanAll(regionStartAbs==binsStart(i) & regionSpanAll==binsRegionSpan(j)));
        end
    end
end
% fsBins(round(end/2), end) = 1; 
fsBins(1, end) = 1;% before ablation: start 0, span all, full pattern
figure, imagesc(binsRegionSpan, binsStart, fsBins, 'alphadata', ~isnan(fsBins));
ylabel('regionStart')
xlabel('regionSpan (seg)')
set(gca,'XAxisLocation','top');
colorbar
axis equal
xlim([min(binsRegionSpan)-0.5, max(binsRegionSpan)+0.5]);
ylim([min(binsStart)-0.5, max(binsStart)+0.5]);
print( [plotDir, 'AblationLocSize'], '-dpdf', '-r0');

tagVal = ~isnan(fsBins);
[xx, yy] = meshgrid(binsRegionSpan, binsStart);
fsBinsInt = griddata(xx(tagVal), yy(tagVal), fsBins(tagVal), xx, yy);

figure, imagesc(binsRegionSpan, binsStart, fsBinsInt, 'alphadata', ~isnan(fsBinsInt))
ylabel('regionStart')
xlabel('regionSpan (seg)')
set(gca,'XAxisLocation','top');
colorbar
axis equal
xlim([min(binsRegionSpan)-0.5, max(binsRegionSpan)+0.5]);
ylim([min(binsStart)-0.5, max(binsStart)+0.5]);
print( [plotDir, 'AblationLocSize_Interpolated'], '-dpdf', '-r0');