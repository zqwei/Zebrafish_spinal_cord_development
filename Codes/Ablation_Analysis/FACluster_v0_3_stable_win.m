%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.  Covariance analysis: Factor analysis -- number of island using non
% LONO methods with active neurons
%
%
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
%


function FACluster_v0_3_stable_win(nFile, thresLMat,thresPsi,thresLifeTime)
addpath('../Func');
setDir;
fileName          = fileNames{nFile}; %#ok<USENS>

load([tempDatDir, fileName, '.mat'], 'new_y');
load([tempDatDir, 'LONOLoading_' fileName, '.mat'], 'CorrectedLMat', 'PsiMat');

side = new_y>0;
[factorComp, factorSizes] =  getFactorIndex(CorrectedLMat, PsiMat, side, thresLMat,thresPsi,thresLifeTime);
save([tempDatDir, 'LONOLoading_' fileName, '_v2.mat'], 'factorComp', 'factorSizes');
end

function [factorComposition, factorSizes] = getFactorIndex(CorrectedLMat, PsiMat, side, thresLMat,thresPsi,thresLifeTime)

factorComposition = cell(numel(CorrectedLMat), 1);
factorSizes       = zeros(numel(CorrectedLMat), 1);
preLMat           = nan(size(CorrectedLMat{1}, 1), 1);

for period = 1:numel(CorrectedLMat)
    LMat          = CorrectedLMat{period};
    currentPsiMat = PsiMat{period};
    LMat(isnan(LMat)) = 0;
    LMat(:, sum(LMat, 1)==0) = [];
    
    % threshold LMat based on lambda and psi threshold
    LMat(currentPsiMat>thresPsi, :) = 0;
    LMat(abs(LMat)<thresLMat) = 0;
    
    % remove bilateral factors
    sideCount = zeros(size(LMat, 2), 1);
    for nFactor = 1:size(LMat, 2)
        sideCount(nFactor) = numel(unique(side(LMat(:, nFactor)>0)));
    end
    LMat(:, sideCount>1) = [];
    
    % remove the single-unit factor in movie
    LMat(:, sum(LMat>0, 1)<=1) = []; % drop factors with zero weight
    
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
            end
        end
    end
    for nFactor     = 1:size(LMat, 2)
        factorComposition{period, factorIndex(nFactor)} = find(LMat(: , nFactor)>0);
        factorSizes(period, factorIndex(nFactor)) = sum(LMat(: , nFactor)>0);
    end
end
% get rid of short-lived factors - those only exist for 1 time window
validFactor = sum(factorSizes>0, 1)>thresLifeTime;
factorComposition = factorComposition(:, validFactor);
factorSizes = factorSizes(:, validFactor);
end