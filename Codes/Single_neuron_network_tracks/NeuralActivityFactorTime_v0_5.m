%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LMatCorrection_v_0_1
% PreLMatTracker_v_0_1
% Spectrogram_v0_3
% Spectrogram_v0_4
%
% Similar to NeuralActivityFactorTime_v0_2 but using new_activeNeuronMat
%
% -------------------------------------------------------------------------
% 
% Ziqiang Wei
% weiz@janelia.hhmi.org
%
%

function NeuralActivityFactorTime_v0_5(nFile)
    addpath('../Func');
    setDir;    
    fileName   = fileNames{nFile};   %#ok<*USENS>
    load([tempDatDir, fileName, '.mat'], 'dff', 'timePoints', 'new_x', 'new_y', 'new_z')
    load([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_1.mat'], 'new_activeNeuronMat', 'oscNeuronMat') 
    load([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_2.mat'], 'preLMat', 'preLMatIndex', 'preLMatTime') 
    numTime           = size(new_activeNeuronMat, 2);
    numNeuron         = size(new_activeNeuronMat, 1);
    neuronFactorIndex = false(size(new_activeNeuronMat));
    timeWin           = 40;
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
    neuronActTime = nan(numNeuron, 1);
    neuronOscTime = nan(numNeuron, 1);
    
    for nNeuron = 1:numNeuron
        if sum(neuronTimeValue(nNeuron, :)) > 0
            neuronTime(nNeuron) = find(neuronTimeValue(nNeuron, :), 1, 'first');
            time_Start = max(1, neuronTime(nNeuron)-timeWin);
            neuronActTime(nNeuron) = mean(new_activeNeuronMat(nNeuron, time_Start:neuronTime(nNeuron)));
            neuronOscTime(nNeuron) = mean(oscNeuronMat(nNeuron, time_Start:neuronTime(nNeuron)));
%             if neuronOscTime(nNeuron)<=0
%                 neuronOscTime(nNeuron) = nan;
%             end
        end
    end
    
    % activation level distribution
    figure;
    [fout, xout] = ksdensity(neuronActTime(~isnan(neuronActTime)), 0:0.02:1);
    plot(xout, fout, 'linewid', 2)
    xlim([0 1])
    box off
    setPrint(8, 6, [plotNetDir 'SingleNeuronActLevelFactorTime_' fileName])
    
    % oscillation level distribution
    figure;
    [fout, xout] = ksdensity(neuronOscTime(~isnan(neuronOscTime)), 0:0.02:1);
    plot(xout, fout, 'linewid', 2)
    xlim([0 1])
    box off
    setPrint(8, 6, [plotNetDir 'SingleNeuronOscLevelFactorTime_' fileName])
    
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
    
    % scatter plot activation level
    figure;
    hold on;
    scatter(new_x, new_y, neuronActTime*100+0.1, neuronTime, 'filled')
    xlim([0 ceil(max(new_x))])
    ylim([-2 2])
    gridxy(1:ceil(max(new_x)), 0, 'color', 'k', 'linestyle', '--')
    box off
    setPrint(8, 6, [plotNetDir 'SingleNeuronFactorTime_' fileName])    
    
    % scatter plot oscillation level
    figure;
    hold on;
    scatter(new_x, new_y, neuronOscTime*100+0.1, neuronTime, 'filled')
    xlim([0 ceil(max(new_x))])
    ylim([-2 2])
    gridxy(1:ceil(max(new_x)), 0, 'color', 'k', 'linestyle', '--')
    box off
    setPrint(8, 6, [plotNetDir 'SingleNeuronFactorTimeOsc_' fileName])
    
% % %     figure;
% % %     plot(neuronActTime, neuronOscTime, 'o')

    save([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_1.mat'], 'neuronTime', 'neuronActTime', 'neuronOscTime', '-append')
    
    
    % plot of neuron in local communities
    numFactor     = max(preLMatIndex);
    preLMatNeuron = sum(preLMat);
    nStep         = 10;
    existLMat     = false(numNeuron, 1);
    mColor        = cbrewer('qual', 'Dark2',  numFactor, 'cubic');

    figure;
    totNeuron     = 0;
    for nFactor   = 1:numFactor
        timeInd   = preLMatTime(:, preLMatNeuron<4 & preLMatIndex == nFactor);
        zeroTime  = min(timeInd);
        endTime   = min(max(timeInd), zeroTime+timeWin);
        LMat      = preLMat(:, preLMatNeuron<4 & preLMatIndex == nFactor & preLMatTime<=endTime);
        LMatInd   = sum(LMat, 2)>0;
        LMatInd(existLMat & LMatInd) = false;
        existLMat = existLMat | LMatInd;
        zeroTime  = min(timeInd);
        endTime   = min(max(timeInd), zeroTime+timeWin);
        minTime   = max(zeroTime - timeWin, 1);
        timeRange = timePoints(minTime)+1:timePoints(endTime);
        timeMarks = timeRange - timePoints(zeroTime);
        if sum(LMatInd) >0 
            dffValue  = bsxfun(@plus, zscore(dff(LMatInd, timeRange), [], 2), (1:sum(LMatInd))'*nStep + totNeuron*nStep);
            dffMedian = median(dffValue, 2);
            hold on
            plot(timeMarks, dffValue, 'Color', mColor(nFactor,:), 'linewidth', 0.5);
            LMatInd = find(LMatInd);
            for mNeuron = 1:length(LMatInd)
                nNeuron = LMatInd(mNeuron);
                text(timePoints(timeWin)+1, mNeuron*nStep + totNeuron*nStep, [num2str(nNeuron) ';' num2str(neuronOscTime(nNeuron)>0.9) '; ' num2str(neuronActTime(nNeuron),'%.2f')])
            end
            xlim([-timePoints(timeWin) timePoints(timeWin)+10])
            gridxy([], dffMedian+1.5, 'color', 'k', 'linestyle', '--')
            axis off
            totNeuron = totNeuron + length(LMatInd);
        end
    end    
    ylim([-1 (totNeuron+1)*nStep+1])
    gridxy(0, [], 'color', 'k', 'linestyle', '--')
    setPrint(20, 4*numFactor, [plotNetDir 'SingleNeuronDynamicsLocalCommunity_' fileName], 'pdf')
    
    
    % plot all sorted neuron
    [neuronTime, neuronTimeInd] = sort(neuronTime, 'ascend');
    neuronFactor = neuronFactor(neuronTimeInd);
    neuronActTime = neuronActTime(neuronTimeInd)>0.6;
    neuronOscTime = neuronOscTime(neuronTimeInd)>0.6;

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
            text(timePoints(timeWin)+1, nNeuron*4, [num2str(nNeuronTime) '; ' num2str(neuronFactor(nNeuron)) '; ' num2str(neuronOscTime(nNeuron)) '; ' num2str(neuronActTime(nNeuron))])
        end
    end
    
    ylim([-1 sum(~isnan(neuronTime))*4+1])
    gridxy(0,[], 'color', 'k', 'linestyle', '--')
    xlim([-timePoints(timeWin) timePoints(timeWin)+10])
    box off
    xlabel('Time')
    ylabel('dff')
%     title(fileName)
    
    setPrint(20, 40, [plotNetDir 'SingleNeuronDynamics_' fileName], 'tiff')
end