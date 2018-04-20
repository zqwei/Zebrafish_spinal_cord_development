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

function NeuralActivityFactorTime_v0_5_1(nFile)
    addpath('../Func');
    setDir;    
    fileName   = fileNames{nFile};   %#ok<*USENS>
    load([tempDatDir, fileName, '.mat'], 'dff', 'timePoints', 'new_x', 'new_y', 'new_z', 'activeNeuronMat')
    load([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_1.mat'], 'new_activeNeuronMat', 'oscNeuronMat') 
    load([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_2.mat'], 'preLMat', 'preLMatIndex', 'preLMatTime') 
    numTime           = size(new_activeNeuronMat, 2);
    numNeuron         = size(new_activeNeuronMat, 1);
    neuronFactorIndex = false(size(new_activeNeuronMat));
    timeBin           = 31;

    for nTime = 1:numTime
        LMat  = preLMat(:, preLMatTime==nTime);
        if ~isempty(LMat)
            neuronFactorIndex(sum(LMat, 2)>0, nTime) = true;
        end
    end

    neuronTimeValue = neuronFactorIndex;
%     neuronTimeValue = double(neuronFactorIndex & new_activeNeuronMat);
%     for nNeuron     = 1:numNeuron
%         smooth_neuronTimeValue      = smooth([zeros(1, 11), neuronTimeValue(nNeuron, :)], 21);
%         neuronTimeValue(nNeuron, :) = smooth_neuronTimeValue(1:end-11) > 0.3;
%     end
    
%     neuronTimeValue = neuronFactorIndex & new_activeNeuronMat & neuronTimeValue;
    
    
    neuronTime    = nan(numNeuron, 1);
    neuronActTime = nan(numNeuron, 1);
    neuronFASize  = nan(numNeuron, 1);
    
    for nNeuron = 1:numNeuron
        if sum(neuronTimeValue(nNeuron, :)) > 0
            neuronTime(nNeuron) = find(neuronTimeValue(nNeuron, :), 1, 'first');
            actCurrNeuron       = activeNeuronMat(nNeuron, :); %#ok<*NODEF>
            actCurrNeuron       = smooth(double(actCurrNeuron), timeBin);
%             % dataset 15 FA not accurate, correct for a single cell 
%             if nNeuron == 26 && nFile == 15
%                 neuronTime(nNeuron) = 99;
%             end
            if neuronTime(nNeuron)-(timeBin+1)/2 > 0
                neuronActTime(nNeuron) = max(actCurrNeuron(1:neuronTime(nNeuron)-(timeBin+1)/2));
                LMat                   = preLMat(:, preLMatTime==neuronTime(nNeuron));
                if ~(nNeuron == 26 && nFile == 15)
                    neuronFASize(nNeuron)  = max(sum(LMat(:, LMat(nNeuron, :)==1)));
                end
            end
        end
    end
    
    neuronActTimeThres = 0.9;
    
    % activation level distribution
    figure;
    subplot(1, 3, 1)
    [fout, xout] = hist(neuronActTime(~isnan(neuronActTime)), 0:0.05:1);
    stairs(xout, fout, 'linewid', 2)
    xlim([0 1.1])
    box off
    
    subplot(1, 3, 2)
    [fout, xout] = hist(neuronFASize(neuronActTime>neuronActTimeThres), 1:nanmax(neuronFASize));
    stairs(xout, fout, 'linewid', 2)
%     xlim([0 1.1])
    box off
    
    subplot(1, 3, 3)
    [fout, xout] = hist(neuronFASize(neuronActTime<neuronActTimeThres), 1:nanmax(neuronFASize));
    stairs(xout, fout, 'linewid', 2)
%     xlim([0 1.1])
    box off
    
    setPrint(8*3, 6, [plotNetDir 'SingleNeuronActLevelFactorTime_' fileName],'pdf')
    
    save([tempDatNetDir, 'NeuronActTime', fileName, '.mat'], 'neuronActTime', 'neuronTime', 'neuronFASize');

end