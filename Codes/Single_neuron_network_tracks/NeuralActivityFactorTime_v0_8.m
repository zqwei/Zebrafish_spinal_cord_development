%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LMatCorrection_v_0_1
% PreLMatTracker_v_0_1
% Spectrogram_v0_3
% Spectrogram_v0_4
%
% Following NeuralActivityFactorTime_v0_5 with labeled neuron
% 1. local community pacemaker
% 2. local community non-active neuron
% 3. recruiting cells
%
% -------------------------------------------------------------------------
% 
% Ziqiang Wei
% weiz@janelia.hhmi.org
%
%


addpath('../Func');
setDir;  

% LocalNeuronIndSet = [];
% neuronActLevelSet = [];

for nFile = [3 7 10 11 16]
    fileName   = fileNames{nFile};   %#ok<*USENS>
    load([tempDatDir, fileName, '.mat'], 'dff', 'timePoints', 'new_x', 'new_y', 'new_z')
    load([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_1.mat'], 'new_activeNeuronMat', 'oscNeuronMat', 'neuronTime', 'neuronActTime', 'neuronOscTime')
    load([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_2.mat'], 'preLMat', 'preLMatIndex', 'preLMatTime')
    
    numTime           = size(new_activeNeuronMat, 2);
    numNeuron         = size(new_activeNeuronMat, 1);
    neuronFactorIndex = false(size(new_activeNeuronMat));
    timeWin           = 30;
    
    % define neuron identity
    % 1. local community pacemaker
    % 2. local community non-active neuron
    % 3. recruiting cells
    
    mColor = [       0    0.4470    0.7410
            0.9290    0.6940    0.1250
            0.8500    0.3250    0.0980
            0.4940    0.1840    0.5560
            0.4660    0.6740    0.1880
            0.3010    0.7450    0.9330
            0.6350    0.0780    0.1840
            0.4588    0.4392    0.7020
            0.4000    0.4000    0.4000
            0.9059    0.1608    0.5412]; % 2 and 5 are the finally two color
    
    new_activeNeuronMat = double(new_activeNeuronMat);
    
    neuronActLevel      = zeros(numNeuron, 1);
    neuronFactorTime    = nan(numNeuron, 1);
    neuronFactorInd     = zeros(numNeuron, 1);
    
    for nNeuron                   = 1:numNeuron
        neuronTime                = preLMatTime(find(preLMat(nNeuron, :), 1, 'first'));
        if ~isempty(neuronTime)
            neuronFactorInd(nNeuron)  = preLMatIndex(find(preLMat(nNeuron, :), 1, 'first'));
            neuronFactorTime(nNeuron) = neuronTime;
            timeRange                 = neuronTime - (5:min(timeWin, neuronTime-1));
            neuronActLevel(nNeuron)   = mean(new_activeNeuronMat(nNeuron, timeRange));
        end
    end
    
    numLMat = max(preLMatIndex);
    LMats    = false(numNeuron, numLMat);

    for nLMat = 1:numLMat
        nTime = find(preLMatIndex == nLMat, 1, 'first');
        LMats(:, nLMat) = preLMat(:, nTime);
        neuronFactorInd(preLMat(:, nTime)) = nLMat;
        neuronFactorTime(preLMat(:, nTime))= preLMatTime(nTime);
    end
    
    LocalNeuronInd = sum(LMats, 2) > 0;
    
    save([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_3.mat'], 'LocalNeuronInd', 'neuronFactorInd', 'neuronFactorTime', 'LMats', 'neuronActLevel')
    
%     LocalNeuronIndSet = [LocalNeuronIndSet; LocalNeuronInd];
%     neuronActLevelSet = [neuronActLevelSet; neuronActLevel];
    
%     figure;
%     hold on
%     for nNeuron = 1:numNeuron
%         if ~isnan(neuronFactorTime(nNeuron))
%             if LocalNeuronInd(nNeuron)
%                 plot(neuronFactorTime(nNeuron)/60, neuronActLevel(nNeuron), 'ow', 'MarkerFaceColor', mColor(neuronFactorInd(nNeuron), :))
%             else
%                 plot(neuronFactorTime(nNeuron)/60, neuronActLevel(nNeuron), 'sw', 'MarkerFaceColor', mColor(neuronFactorInd(nNeuron), :))
%             end
%         end
%     end
%     hold off
%     xlabel('Time (hour)')
%     ylabel('Active Level')
%     set(gca, 'TickDir', 'out')
%     setPrint(8, 6, [plotNetDir 'SingleNeuronActLevelFactorTimeCellType_' fileName], 'pdf')
%          
%     figure;
%     hold on
%     for nType = 1:2
%         [fout, xout] = ksdensity(neuronActLevel(LocalNeuronInd == nType-1), 0:0.02:1.02, 'Bandwidth', 0.01);
%         stairs(xout, fout/sum(fout), 'linewid', 2)
%         xlim([0 1.01])
%     end
%     box off
%     xlabel('Active level')
%     ylabel('Frac.')
%     setPrint(8, 6, [plotNetDir 'SingleNeuronActLevelFactorTimeCellType_' fileName])

end


% across animal summary
% figure;
% hold on
% for nType = 1:2
%     [fout, xout] = ksdensity(neuronActLevelSet(LocalNeuronIndSet == nType-1), 0:0.02:1.02, 'Bandwidth', 0.01);
%     stairs(xout, fout/sum(fout), 'linewid', 2)
%     xlim([0 1.01])
% end
% box off
% xlabel('Active level')
% ylabel('Frac.')
% setPrint(8, 6, [plotNetDir 'SingleNeuronActLevelFactorTimeCellType_' fileName])