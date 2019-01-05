% 0. evaluate the effect of different threshold on all datasets
%
% This code is for testing only
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%
addpath('../Func');
setDir;

p = 0.05;
nFileList = 25:2:58;

fishListCutA = [1, 2, 3, 4, 7]; % anterior cut
fishListCutM = [8, 9, 10, 12]; % middle cut
fishListCutP = [13, 15, 16, 17]; % posterior cut
fishList = [fishListCutA, fishListCutM, fishListCutP]; % anterior cut
actThresList = 0:0.1:1;
percOverlap = zeros(numel(fishList), numel(actThresList));
for i = 1:numel(fishList)
    
    nFish = fishList(i);
    
    for nExp = 1:2
        nFile = nFileList(nFish) + nExp - 1;
        fileDirName   = fileDirNames{nFile}; %#ok<*USENS>
        fileName      = fileNames{nFile};
        
        dirImageData  = [fileDirName '/'];
        load([dirImageData, 'profile.mat'], 'segAblation');
        if nExp == 1
            load([tempDatDir, fileName, '.mat'], 'activeNeuronMat');
            activeTagBefore = activeNeuronMat;
        else
            load([tempDatDir, fileName, '_full.mat'], 'activeNeuronMat');
        end
    end
    
    for j = 1:numel(actThresList)
        activeThres = actThresList(j);
        activeTag = sum(activeNeuronMat, 2)>= 10* activeThres;
%         % maximize the overlap of active neuron population - not meaningful
%         percOverlap(i, j) = sum(activeTag & activeTagBefore)/sum(activeTag | activeTagBefore); 
%         % maximize the number of neurons with same activation-classification
%         percOverlap(i, j) = sum(~xor(activeTag, activeTagBefore))/numel(activeTag); 
        % find the cutoff that best distinguish the active-before and nonactive-before set
        percOverlap(i, j) = 2 - sum(~activeTag & activeTagBefore)/sum(activeTagBefore) - sum(activeTag & ~activeTagBefore)/sum(~activeTagBefore) ; 
    end
end

figure, imagesc(percOverlap);
[maxOverlap, thresInd] = max(percOverlap,[],  2);
hold on;
text(thresInd, 1:numel(fishList), '*');
useThres = actThresList(thresInd);
colorbar;
