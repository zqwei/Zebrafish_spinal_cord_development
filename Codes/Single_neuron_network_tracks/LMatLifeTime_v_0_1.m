%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LMatCorrection_v_0_1
% 
% based on LMatCorrection_v_0_1 corrected L-Matrix:
% compute the life time of each factor with in timeWin = 20 min
% aim to remove the ones with short life time
%
% -------------------------------------------------------------------------
% 
% Ziqiang Wei
% weiz@janelia.hhmi.org
%
%


function LMatLifeTime_v_0_1(nFile)            
    addpath('../Func');
    setDir;    
    fileName   = fileNames{nFile};   %#ok<*USENS>
    load([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_1.mat'], 'CorrectedLMat') 
    numTime    = length(CorrectedLMat);
    timeWin    = 20;
    lifeTimeTable = CorrectedLMat;
    totFactor  = 0;

    for nTime = 1:numTime - 1
        LMat  = CorrectedLMat{nTime};
        if ~isempty(LMat)
            LMat             = LMat > 0; % make LMat as a index only matrix; point set operations
            lifeTime         = zeros(1, size(LMat, 2));
            timeLen          = min(timeWin, numTime-nTime);
            for mTime        = (nTime + 1): (nTime + timeLen)
                LMatTime     = CorrectedLMat{mTime};
                if ~isempty(LMatTime)
                    LMatTime        = LMatTime > 0;
                    similarityScore = max(double(LMat)' * double(LMatTime), [], 2) ./ sum(LMat)';
                    lifeTime(similarityScore > 0.4999) = lifeTime(similarityScore > 0.4999) + 1;
                end
            end
            lifeTimeTable{nTime} = lifeTime / timeLen;
        else
            lifeTimeTable{nTime} = 0;
        end
    end
    
    LMat                   = CorrectedLMat{numTime};
    lifeTimeTable{numTime} = ones(1, size(LMat, 2));
    save([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_1.mat'], 'lifeTimeTable', '-append') 
end