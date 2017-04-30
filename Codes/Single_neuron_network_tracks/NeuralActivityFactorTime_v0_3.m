%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% based on result in NeuralActivityFactorTime_v0_2
% computation from PreLMatTracker_v_0_1
%
% -------------------------------------------------------------------------
% 
% Ziqiang Wei
% weiz@janelia.hhmi.org
%
%

function NeuralActivityFactorTime_v0_3(nFile)
    addpath('../Func');
    setDir;    
    fileName   = fileNames{nFile};   %#ok<*USENS>
    load([tempDatDir, fileName, '.mat'], 'dff', 'timePoints', 'new_x', 'new_y', 'new_z')
    load([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_1.mat'], 'new_activeNeuronMat') 
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
    
    for nNeuron = 1:numNeuron
        if sum(neuronTimeValue(nNeuron, :)) > 0
            neuronTime(nNeuron) = find(neuronTimeValue(nNeuron, :), 1, 'first');
            neuronActTime(nNeuron) = mean(new_activeNeuronMat(nNeuron, max(1, (neuronTime(nNeuron)-timeWin):neuronTime(nNeuron))));
        end
    end
    
    figure;
    [fout, xout] = ksdensity(neuronActTime(~isnan(neuronActTime)), 0:0.02:1);
    plot(xout, fout, 'linewid', 2)
    xlim([0 1])
    box off
    setPrint(8, 6, [plotNetDir 'SingleNeuronActLevelFactorTime_' fileName])
    
end