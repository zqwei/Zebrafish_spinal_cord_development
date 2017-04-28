%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LMatCorrection_v_0_1
% PreLMatTracker_v_0_1
% 
% based on LMatCorrection_v_0_1 corrected L-Matrix and lifetime factor
% computation from PreLMatTracker_v_0_1
%
% -------------------------------------------------------------------------
% 
% Ziqiang Wei
% weiz@janelia.hhmi.org
%
%

function NeuralActivityFactorTime_v0_2(nFile)
    addpath('../Func');
    setDir;    
    fileName   = fileNames{nFile};   %#ok<*USENS>
    load([tempDatDir, fileName, '.mat'], 'dff', 'timePoints', 'new_x', 'new_y', 'new_z')
    load([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_1.mat'], 'new_activeNeuronMat') 
    load([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_2.mat'], 'preLMat', 'preLMatIndex', 'preLMatTime') 
    numTime           = size(new_activeNeuronMat, 2);
    numNeuron         = size(new_activeNeuronMat, 1);
    neuronFactorIndex = false(size(new_activeNeuronMat));
    timeWin           = 30;
    linew             = 0.5;
    colorSet          = cbrewer('qual', 'Dark2',  numTime, 'cubic');

    for nTime = 1:numTime
        LMat  = preLMat(:, preLMatTime==nTime);
        if ~isempty(LMat)
            neuronFactorIndex(sum(LMat, 2)>0, nTime) = true;
        end
    end
    
    neuronTimeValue = double(neuronFactorIndex & new_activeNeuronMat);
    for nNeuron     = 1:numNeuron
        smooth_neuronTimeValue      = smooth([zeros(1, 11), neuronTimeValue(nNeuron, :)], 21);
        neuronTimeValue(nNeuron, :) = smooth_neuronTimeValue(1:end-11) > 0.3;
    end
    
    neuronTimeValue = neuronFactorIndex & new_activeNeuronMat & neuronTimeValue;
    
    
    neuronTime = nan(numNeuron, 1);
    
    for nNeuron = 1:numNeuron
        if sum(neuronTimeValue(nNeuron, :)) > 0
            neuronTime(nNeuron) = find(neuronTimeValue(nNeuron, :), 1, 'first');
        end
    end
    
    neuronFactor = nan(numNeuron, 1);
    for nNeuron = 1:numNeuron
        if ~isnan(neuronTime(nNeuron))
            LMat                  = preLMat(:, preLMatTime == neuronTime(nNeuron));
            LMatIndex             = preLMatIndex(:, preLMatTime == neuronTime(nNeuron));
            LMatNeuron            = LMat(nNeuron, :) >0;
            LMat                  = LMat(:, LMatNeuron);
            LMatIndex             = LMatIndex(:, LMatNeuron);
            [~, largeFactor]      = max(sum(LMat, 1));
            if length(largeFactor) ~= 1; keyboard(); end
            neuronFactor(nNeuron) = LMatIndex(largeFactor);
        end
    end
    
    figure;
    hold on;
    scatter(new_x, new_y, [], neuronTime, 'filled')
    xlim([0 ceil(max(new_x))])
    ylim([-2 2])
    gridxy([0:ceil(max(new_x))], [0], 'color', 'k', 'linestyle', '--')
    box off
    title(fileName)

    [neuronTime, neuronTimeInd] = sort(neuronTime, 'ascend');
    neuronFactor = neuronFactor(neuronTimeInd);

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
            text(timePoints(timeWin)+1, nNeuron*4, [num2str(nNeuronTime) ';' num2str(neuronFactor(nNeuron))])
        end
    end
    
    ylim([-1 sum(~isnan(neuronTime))*4+1])
    gridxy([0],[], 'color', 'k', 'linestyle', '--')
    xlim([-timePoints(timeWin) timePoints(timeWin)+10])
    box off
    xlabel('Time')
    ylabel('dff')
    title(fileName)
%     
%     setPrint(20, 40, [plotNetDir 'SingleNeuronDynamics_' fileName], 'tiff')
end