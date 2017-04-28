%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LMatCorrection_v_0_1
% LMatLifeTime_v_0_1
% 
% based on LMatCorrection_v_0_1 corrected L-Matrix and lifetime factor
% computation froom LMatLifeTime_v_0_1
%
% -------------------------------------------------------------------------
% 
% Ziqiang Wei
% weiz@janelia.hhmi.org
%
%

function NeuralActivityFactorTime_v0_1(nFile)
    addpath('../Func');
    setDir;    
    fileName   = fileNames{nFile};   %#ok<*USENS>
    load([tempDatDir, fileName, '.mat'], 'dff', 'timePoints')
    load([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_1.mat'], 'CorrectedLMat', 'lifeTimeTable', 'new_activeNeuronMat') 
    numTime           = length(CorrectedLMat);
    numNeuron         = size(new_activeNeuronMat, 1);
    neuronFactorIndex = false(size(new_activeNeuronMat));
    lifeTimeThres     = 0.5;
    timeWin           = 20;
    linew             = 0.5;
    colorSet          = cbrewer('seq', 'PuBuGn',  numTime, 'cubic');

    for nTime = 1:numTime
        LMat  = CorrectedLMat{nTime};
        if ~isempty(LMat)
            lifeTime = lifeTimeTable{nTime};
            LMat     = LMat(:, lifeTime > lifeTimeThres);
            neuronFactorIndex(sum(LMat, 2)>0, nTime) = true;
        end
    end
    
    neuronTimeValue = neuronFactorIndex & new_activeNeuronMat;
    neuronTime = nan(numNeuron, 1);
    
    for nNeuron = 1:numNeuron
        if sum(neuronTimeValue(nNeuron, :)) > 0
            neuronTime(nNeuron) = find(neuronTimeValue(nNeuron, :), 1, 'first');
        end
    end
    
    [neuronTime, neuronTimeInd] = sort(neuronTime, 'ascend');
    figure;
    hold on
    for nNeuron = 1:length(neuronTime)
        nNeuronTime = neuronTime(nNeuron);
        nNeuronInd  = neuronTimeInd(nNeuron);
        
        if ~isnan(nNeuronTime)    
            minTime   = max(nNeuronTime - timeWin, 1);
            maxTime   = min(numTime, nNeuronTime + timeWin);
            timeRange = timePoints(minTime)+1:timePoints(maxTime);
            timeMarks = timeRange - timePoints(nNeuronTime);
            dffValue  = zscore(dff(nNeuronInd, timeRange))+ nNeuron*4;
            plot(timeMarks, dffValue, 'Color', colorSet(nNeuronTime,:), 'linewidth', linew)
        end
    end
    
    ylim([-1 sum(~isnan(neuronTime))*4+1])
%     ylim([-1 10*4+1])
    xlim([-timePoints(timeWin) timePoints(timeWin)])
    box off
    xlabel('Time')
    ylabel('dff')
    
    setPrint(20, 40, [plotNetDir 'SingleNeuronDynamics_' fileName], 'tiff')
end