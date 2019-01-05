%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.1 (Single cut) FA after cut, all fish, matrix view
%
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%

function Ablation_v_3_1(fishListCutA, fishListCutM, fishListCutP)
addpath('../Func');
setDir;



fishList = [fishListCutA, fishListCutM, fishListCutP];

nFileList = 25:2:100;

fracPatternNeuron = zeros(numel(fishList), 4);
factorSizeAll = nan(numel(fishList)*2, 2);
factorSpanAll = nan(numel(fishList)*2, 2);

fsSpanThresFull = 0.5;
fsSpanThresFullPart = 0.05;

for i = 1:numel(fishList)
    nFish = fishList(i);
    nFile = nFileList(nFish) + 1;
    fileDirName   = fileDirNames{nFile}; %#ok<*USENS>
    fileName      = fileNames{nFile};
    
    dirImageData  = [fileDirName '/'];
    load([dirImageData, 'profile.mat'], 'segAblation');
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
        nLR = mod(nSec, 2);
        nAP = floor((nSec+1)/2);
        factorSizeAll(i*2+nLR-1, nAP) = max(factorSizes);
        factorSpanAll(i*2+nLR-1, nAP) = max(factorSpans);
    end
    disp('');
end


m = factorSpanAll;
m(m>fsSpanThresFull) = 2;
m(m<fsSpanThresFull & m>fsSpanThresFullPart) = 1;
m(m<fsSpanThresFullPart) = 0;
im = plotMatrixType(m', 1);
imwrite(im, [plotDir, 'SingleAblationMatrix_HalfFish.tif']);

% m = (factorSpanAll(1:2:end, :) + factorSpanAll(2:2:end, :))/2;
m = max(factorSpanAll(1:2:end, :), factorSpanAll(2:2:end, :));
m(m>fsSpanThresFull) = 2;
m(m<fsSpanThresFull & m>fsSpanThresFullPart) = 1;
m(m<fsSpanThresFullPart) = 0;
im = plotMatrixType(m', 0);
imwrite(im, [plotDir, 'SingleAblationMatrix_WholeFish.tif']);

