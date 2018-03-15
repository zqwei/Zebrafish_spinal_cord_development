%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% PreLMat index
% 
% Track the history of precesor of each factor
% 
% based on code -- 
% Yinan's version of FACluster_v0_5
% 
% -------------------------------------------------------------------------
% 
% Ziqiang Wei
% weiz@janelia.hhmi.org
%
%


function PreLMatTracker_v_0_1(nFile)            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile};   %#ok<*USENS>
    load([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_1.mat'], 'CorrectedLMat', 'new_activeNeuronMat', 'lifeTimeTable') 

    numTime           = length(CorrectedLMat);
    numNeuron         = size(new_activeNeuronMat, 1);
    halfActTime       = nan(numNeuron, 1);
    timeBin           = 15;
    activeThres       = 0.65;
    lifeTimeThres     = 0.5;
    LMatThres         = 0.5;
    timeWin           = 10;
    preLMat           = [];
    preLMatIndex      = [];
    preLMatTime       = [];
    

    for nTime = 1:numTime
        LMat      = CorrectedLMat{nTime};
        activeTag = new_activeNeuronMat(:, nTime);
        lifeTime  = lifeTimeTable{nTime};
        
        % remove the factors with short life time
        LMat      = LMat(:, lifeTime >= lifeTimeThres);
        
        % remove the factors with weak link (in order to remove the overlapped factors)
        LMat      = LMat > LMatThres;
        LMat(:,sum(LMat, 1)<2) = [];
        
        % determine the factor index -- color index part
        if size(LMat,2) >= 1
            if isempty(preLMat)
%                 factorIndex  = 1:size(LMat, 2);
                preLMat      = LMat;
                preLMatIndex = 1:size(LMat, 2);
                totInd       = size(LMat, 2);
                preLMatTime  = ones(1, totInd)*nTime;
            else                
                sizeLMat        = sum(LMat, 1);
                [~, indexLMat]  = sort(sizeLMat, 'descend');
                LMat            = LMat(:, indexLMat);
                similarityScore = double(preLMat') * double(LMat);
                [similarityValue, ind] =  max(similarityScore, [], 1);
                
                diffInd         = similarityValue < 2;
                
                for nFactor     = find(~diffInd)
                    if sum(similarityScore(:, nFactor) == similarityValue(nFactor)) > 1
                        allIndSet         = find(similarityScore(:, nFactor) == similarityValue(nFactor));
                        allIndSize        = sum(preLMat(:, allIndSet));
                        [~, indLargeSize] = max(allIndSize);
                        ind(nFactor)      = allIndSet(indLargeSize);
                    end
                end
                
                LMatIndex       = preLMatIndex(ind);
                LMatIndex(diffInd) = totInd + (1:sum(diffInd));
                preLMat         = [preLMat, LMat];
                preLMatIndex    = [preLMatIndex, LMatIndex];
                totInd          = totInd + sum(diffInd);
                preLMatTime     = [preLMatTime, ones(1, size(LMat, 2))*nTime];                  
            end
        end
    end
    
    save([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_2.mat'], 'preLMat', 'preLMatIndex', 'preLMatTime');
    
end