%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Track factor identity and analyze the composition of initial factors
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%
function [histPairType, actLevelBeforePattern] = Leader_v4(nFile)
% track factor identity using similarity between LMat
addpath('../Func')
setDir;
fileName = fileNames{nFile};
load([tempDatDir '/' fileName '.mat'], 'mnx', 'activeNeuronMat');
load([tempDatDir '/LONOLoading_' fileName '.mat'], 'CorrectedLMat');
load([tempDatDir, 'Leader_', fileName, '.mat'], 'factorSize', 'activeTime');

[factorComp, factorSizes] = getFactorIndex(CorrectedLMat);
figure, imagesc(factorSizes')
xlabel('Time index');
ylabel('Factor index');
setPrint(8, 6, [plotDir,  'factorSizeIdentityEvol_' fileName], 'pdf');
close

leaderThres = 1.5;
nNeuron = numel(factorSize);
pairsCombination = zeros(nNeuron, nNeuron);
actLevelBeforePattern = zeros(0, 2); % activity level (percentage of active window before patterning, 1st&2nd active in the pair)
histPairType = zeros(3, 1);
for nTime = 1:size(factorSizes, 1)
    for nFactor = 1:size(factorSizes, 2)
        if factorSizes(nTime, nFactor)==2 
            pair = factorComp{nTime, nFactor};
            pairsCombination(pair(1), pair(2)) = pairsCombination(pair(1), pair(2))+1;
            pairsCombination(pair(2), pair(1)) = pairsCombination(pair(2), pair(1))+1;
            if pairsCombination(pair(1), pair(2)) == 1 && activeTime(pair(1))>0 && activeTime(pair(2))>0 
                if factorSize(pair(1))<leaderThres && factorSize(pair(2))<leaderThres
                    histPairType(1) = histPairType(1)+1;
                elseif factorSize(pair(1))>=leaderThres && factorSize(pair(2))>=leaderThres
                    histPairType(3) = histPairType(3)+1;
                else
                    histPairType(2) = histPairType(2)+1;
                end
                if activeTime(pair(1))>activeTime(pair(2))
                    pair = pair([2, 1]); % pair index: [fist active neuron, second active neuron]
                end
                firstActiveTime(1) = find(activeNeuronMat(pair(1), :)>0, 1);
                firstActiveTime(2) = find(activeNeuronMat(pair(2), :)>0, 1);

                actLevel = [mean(activeNeuronMat(pair(1), firstActiveTime(1):nTime-1)), mean(activeNeuronMat(pair(2), firstActiveTime(2):nTime-1))];
                actLevel(isnan(actLevel)) = 0;
                actLevelBeforePattern = [actLevelBeforePattern; actLevel];
            end
        end
    end
end


figure, bar(histPairType);
set(gca, 'xticklabel', {'Leader-Leader', 'Leader-NonLeader', 'NonLeader-NonLeader'});
xlabel('Type of Pair')
ylabel('Count');
setPrint(8, 6, [plotDir,  'InitialPairTypes_' fileName], 'pdf');

[histCount, centers] = hist3(actLevelBeforePattern);
figure, imagesc(centers{1}, centers{2}, histCount);
% figure, hist3(actLevelBeforePattern);
set(gca, 'YDir', 'normal');
xlabel('Activity level - 1st in pair')
ylabel('Activity level - 2nd in pair')
zlabel('Number of pairs')
colorbar
[p, ~] = signrank(actLevelBeforePattern(:, 1), actLevelBeforePattern(:, 2), 'tail', 'left');
title(['p=' num2str(p, '%.3f')]);
setPrint(8, 6, [plotDir,  'ActLevelInitialPairs_' fileName], 'pdf');
close;


end

function [factorComposition, factorSizes] = getFactorIndex(CorrectedLMat)
factorComposition = cell(numel(CorrectedLMat), 1);
preLMat           = nan(size(CorrectedLMat{1}, 1), 1);

for period = 1:numel(CorrectedLMat)
    LMat          = CorrectedLMat{period};
    LMat(isnan(LMat)) = 0;
    LMat(:, sum(LMat, 1)==0) = [];
    % drop factors completely contained within another factor
    invIndex = [];
    currLMat = LMat>0;
    for ii = 1:size(LMat, 2)
        for jj = ii+1:size(LMat, 2)
            if sum(ismember(find(currLMat(:, ii)), find(currLMat(:, jj))))==sum(currLMat(:, ii))
                invIndex = [invIndex, ii];
            elseif sum(ismember(find(currLMat(:, ii)), find(currLMat(:, jj))))==sum(currLMat(:, jj))
                invIndex = [invIndex, jj];
            end
        end
    end
    LMat(:, invIndex) = [];
    % determine the factor index
    if size(LMat,2) >= 1
        if sum(~isnan(preLMat(:))) == 0
            factorIndex  = 1:size(LMat, 2);
            preLMat      = LMat;
        else
            
            [~, maxFactorPerNeuronIndex] = max(LMat(sum(LMat, 2)>0, :), [], 2);
            sideRemoveList  = histc(maxFactorPerNeuronIndex, 1:size(LMat, 2)) <2; % remove the factor has no dominate factors
            
            LMat(:, sideRemoveList) = [];
            LMat            = LMat > 0;
            sizeLMat        = sum(LMat, 1);
            [~, indexLMat]  = sort(sizeLMat, 'descend');
            LMat            = LMat(:, indexLMat);
            factorIndex     = zeros(size(LMat, 2), 1);
            
            % compute similarity matrix
            similarityScore = zeros(size(LMat, 2), size(preLMat, 2));
            for nFactor     = 1:size(LMat, 2)
                if sum(isnan(LMat(:)))>0; keyboard();end
                similarityScore(nFactor, :) = sum(bsxfun(@and, LMat(:, nFactor), preLMat));
            end
            
            % check if any prefactor has no connection with new factors
            % decide which factor is not included in preLMatIndex
            % maxIndex is the factor with the maximum coverage with the prefactors
            [~, maxIndex]   = max(similarityScore, [], 1);
            
            % check if any prefactor is merged (factor has maximum coverages with more than one prefactors, pick the larger one as its index)
            for nFactor     = 1:size(LMat, 2)
                nFacotrNumPreFactor = sum(maxIndex == nFactor);
                switch nFacotrNumPreFactor
                    case 0
                        preLMat = [preLMat, LMat(:, nFactor)];
                        factorIndex(nFactor) = size(preLMat, 2);
                    case 1
                        if similarityScore(nFactor, maxIndex == nFactor) == 0
                            preLMat = [preLMat, LMat(:, nFactor)];
                            factorIndex(nFactor) = size(preLMat, 2);
                        else
                            preLMatIndex         = find(maxIndex == nFactor);
                            factorIndex(nFactor) = preLMatIndex;
                            preLMat(:, preLMatIndex) = preLMat(:, preLMatIndex) | LMat(:, nFactor);
                        end
                    otherwise
                        preLMatIndex         = find(maxIndex == nFactor);
                        [~, nFactorMaxIndex] = max(similarityScore(nFactor, preLMatIndex));
                        factorIndex(nFactor) = preLMatIndex(nFactorMaxIndex);
                        preLMat(:, preLMatIndex(nFactorMaxIndex)) = preLMat(:, preLMatIndex(nFactorMaxIndex)) | LMat(:, nFactor);
                end
                factorComposition{period, factorIndex(nFactor)} = find(LMat(: , nFactor)>0);
                factorSizes(period, factorIndex(nFactor)) = sum(LMat(: , nFactor)>0);
            end
        end
    end
end
end