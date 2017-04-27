%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LMatCorrection_v_0_1
% 
% based on EV corrected L-Matrix:
% 1. drop non-active neurons in factors
% 2. drop part of the factor with any side = 1
% 3. set valid threshould of LMat at 0.3
%
%
% -------------------------------------------------------------------------
% 
% Ziqiang Wei
% weiz@janelia.hhmi.org
%
%


function LMatCorrection_v_0_1(nFile)            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile};   %#ok<*USENS>
    load([tempDatDir, fileName, '.mat'], 'side', 'activeNeuronMat', 'timePoints'); 
    load([tempDatDir, 'LONOLoading_' fileName, '.mat'], 'CorrectedLMat') 

    numTime           = length(CorrectedLMat);
    numNeuron         = length(side);
    halfActTime       = nan(numNeuron, 1);
    timeBin           = 15;
    activeThres       = 0.65;
    new_activeNeuronMat = activeNeuronMat;
    LMatThres         = 0.2;
    
    
    % compute half-activation time
    for nNeuron       = 1:numNeuron
        actCurrNeuron   = activeNeuronMat(nNeuron, :); %#ok<*NODEF>
        firstActiveTime = find(smooth(double(actCurrNeuron), timeBin) > activeThres, 1, 'first');
        if ~isempty(firstActiveTime)
            removeIndex     = find(actCurrNeuron==0);
            removeIndex(removeIndex<firstActiveTime) = [];
            removeTimeIndex = false(1, numTime);
            removeTimeIndex(removeIndex) = true; % remove all zeros after first activation time
        else
            removeTimeIndex = false(1, numTime);
        end
        
        if ~isempty(firstActiveTime)
            mdl = fitglm(find(~removeTimeIndex), activeNeuronMat(nNeuron, ~removeTimeIndex), 'linear', 'Distribution', 'binomial');
            beta = mdl.Coefficients.Estimate;
            halfActTime(nNeuron) = -beta(1)/beta(2);
            if (mean(activeNeuronMat(nNeuron, 1:ceil(numTime*.5))) > 0.80 ||  mean(activeNeuronMat(nNeuron, :)) > 0.80) ...
                    && (halfActTime(nNeuron)<0 || halfActTime(nNeuron)>numTime)
                halfActTime(nNeuron) = 0;
            end
            if (halfActTime(nNeuron)<0 || halfActTime(nNeuron)>numTime); halfActTime(nNeuron) = nan; end
        end
        if isnan(halfActTime(nNeuron))
            new_activeNeuronMat(nNeuron, :) = false;
        else
            new_activeNeuronMat(nNeuron, 1:ceil(halfActTime(nNeuron))) = false;
        end
    end

    preLMat           = nan(numNeuron, 1);
    

    for nTime = 1:numel(timePoints)
        LMat          = CorrectedLMat{nTime};
        LMat(isnan(LMat)) = 0;
        LMat(:, sum(LMat, 1)==0) = [];
        activeTag = new_activeNeuronMat(:, nTime);
                
        LMatNeuron    = LMat>LMatThres;
        numFactor     = size(LMatNeuron, 2);
        % drop non-active neurons
        LMatNeuron(~activeTag, :) = false;
        % code to drop part of the factor with any side = 1
        for nFactor   = 1:size(LMatNeuron, 2)
            neuronFactor = LMatNeuron(:, nFactor);
            dominateSide = ceil(mean(side(neuronFactor)) - 0.5);
            otherSide    = 3 - dominateSide;
            if sum(neuronFactor & side == dominateSide) == 1
                LMatNeuron(side == dominateSide, nFactor) = false;
            end
            if sum(neuronFactor & side == otherSide) == 1
                LMatNeuron(side == otherSide, nFactor) = false;
            end
        end
        
        % code to drop overlapped factors
%         for nFactor   = 1:size(LMatNeuron, 2)
%             if sum(LMatNeuron(:,nFactor)) >0 % skip if nFactor is dropped
%                 for mFactor = nFactor+1:size(LMatNeuron, 2)
%                     if sum(LMatNeuron(:,mFactor)) >0 % skip if mFactor is dropped
%                         if all(ismember(find(LMatNeuron(:,nFactor))', find(LMatNeuron(:,mFactor))')) % if nFactor inside any other factors, drop nFactor
%                             LMatNeuron(:,nFactor) = 0;
%                             continue;
%                         elseif all(ismember(find(LMatNeuron(:,mFactor))', find(LMatNeuron(:,nFactor))')) % if mFactor inside any other factors, drop nFactor
%                             LMatNeuron(:,mFactor) = 0;
%                             continue;
%                         end
%                     end
%                 end
%             end
%         end
        % end of drop overlapped factor code
        
        LMat(LMatNeuron == 0)    = 0;
        LMat(:, sum(LMat, 1)==0) = []; % drop factors with zero weight 
        CorrectedLMat{nTime}     = LMat;
    end
    save([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_1.mat'], 'CorrectedLMat', 'new_activeNeuronMat') 
end