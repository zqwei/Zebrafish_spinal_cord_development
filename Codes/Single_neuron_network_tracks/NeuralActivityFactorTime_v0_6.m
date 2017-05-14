%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LMatCorrection_v_0_1
% PreLMatTracker_v_0_1
% Spectrogram_v0_3
% Spectrogram_v0_4
%
% Summary of NeuralActivityFactorTime_v0_5
%
% -------------------------------------------------------------------------
% 
% Ziqiang Wei
% weiz@janelia.hhmi.org
%
%

addpath('../Func');
setDir;   
fileToAnalysis = [3 7 10 16]; % 
mColor         = cbrewer('qual', 'Dark2',  32, 'cubic');
figure;
hold on
barPlot   = nan(32, 4);

for nFile = fileToAnalysis

 
    fileName   = fileNames{nFile};   %#ok<*USENS>
    load([tempDatDir, fileName, '.mat'], 'dff', 'timePoints', 'new_x', 'new_y', 'new_z')
    load([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_1.mat'], 'new_activeNeuronMat', 'oscNeuronMat', 'neuronTime', 'neuronActTime', 'neuronOscTime') 
    load([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_2.mat'], 'preLMat', 'preLMatIndex', 'preLMatTime') 
    numTime           = size(new_activeNeuronMat, 2);
    numNeuron         = size(new_activeNeuronMat, 1);
    neuronFactorIndex = false(size(new_activeNeuronMat));
    timeWin           = 40;


    % plot of neuron in local communities
    numFactor     = max(preLMatIndex);
    preLMatNeuron = sum(preLMat);
    nStep         = 10;
    existLMat     = false(numNeuron, 1);
    LocalMaxSize  = 5;
    timeWinEnd    = 10;

    for nFactor   = 1:numFactor
        timeInd   = preLMatTime(:, preLMatNeuron<LocalMaxSize & preLMatIndex == nFactor);
        if max(timeInd) - min(timeInd) < timeWinEnd; continue; end
        zeroTime  = min(timeInd);
        if isempty(zeroTime); continue; end
        endTime   = min(max(timeInd), zeroTime+timeWinEnd);
        LMat      = preLMat(:, preLMatNeuron<LocalMaxSize & preLMatIndex == nFactor & preLMatTime<=endTime);
        LMatInd   = sum(LMat, 2)>0;
        LMatInd(existLMat & LMatInd) = false;
        existLMat = existLMat | LMatInd;
        numOscNeuron = sum(neuronOscTime(LMatInd)>0.75 & neuronActTime(LMatInd)>0.75);
        if sum(LMatInd)>1
            if isnan(barPlot(nFile, numOscNeuron+1))
                barPlot(nFile, numOscNeuron+1) = 1;
            else
                barPlot(nFile, numOscNeuron+1) = barPlot(nFile, numOscNeuron+1)+1;
            end
        end
    end

end

bar(0:3, barPlot', 'stacked'), colormap(mColor)
box off
xlabel('Number oscillatory neurons')
ylabel('Number local communities')
xlim([-0.5 2.5])
set(gca, 'TickDir', 'Out')
setPrint(8, 6, [plotNetDir 'SingleNeuronDynamicsLocalCommunity_summary'], 'pdf')



