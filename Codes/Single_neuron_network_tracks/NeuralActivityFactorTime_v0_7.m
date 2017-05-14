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

function NeuralActivityFactorTime_v0_7(nFile)
    addpath('../Func');
    setDir;    
    fileName   = fileNames{nFile};   %#ok<*USENS>
    load([tempDatDir, fileName, '.mat'], 'dff', 'timePoints', 'new_x', 'new_y', 'new_z')
    load([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_1.mat'], 'new_activeNeuronMat', 'oscNeuronMat', 'neuronTime', 'neuronActTime', 'neuronOscTime') 
    load([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_2.mat'], 'preLMat', 'preLMatIndex', 'preLMatTime') 
    numTime           = size(new_activeNeuronMat, 2);
    numNeuron         = size(new_activeNeuronMat, 1);
    neuronFactorIndex = false(size(new_activeNeuronMat));
    timeWin           = 40;
    linew             = 0.5;
    colorSet          = cbrewer('qual', 'Dark2',  numTime, 'cubic');
    
    % define neuron identity
    % 1. local community pacemaker
    % 2. local community non-active neuron
    % 3. recruiting cells

    numFactor     = max(preLMatIndex);
    preLMatNeuron = sum(preLMat);
    LocalMaxSize  = 4;
    timeWinEnd    = 10;
    neuronType    = nan(numNeuron, 1);
    neuronType(sum(preLMat, 2)>0) = 3;
    oscInd        = neuronOscTime>0.75 & neuronActTime>0.75;
    
    for nFactor   = 1:numFactor
        timeInd   = preLMatTime(:, preLMatNeuron<LocalMaxSize & preLMatIndex == nFactor);
        if max(timeInd) - min(timeInd) < timeWinEnd; continue; end
        zeroTime  = min(timeInd);
        if isempty(zeroTime); continue; end
        endTime   = min(max(timeInd), zeroTime+timeWinEnd);
        LMat      = preLMat(:, preLMatNeuron<LocalMaxSize & preLMatIndex == nFactor & preLMatTime<=endTime);
        LMatInd   = sum(LMat, 2)>0;
        neuronType(LMatInd & oscInd) = 1;
        neuronType(LMatInd & ~oscInd & neuronType~=1) = 2;
    end    
    
    figure;
    hold on
    for nType = 1:3
        [fout, xout] = ksdensity(neuronActTime(neuronType == nType), 0:0.02:1, 'Bandwidth', 0.02);
        stairs(xout, fout/sum(fout), 'linewid', 2)
        xlim([0 1])
    end
    box off
    xlabel('Active level')
    ylabel('Frac.')
    setPrint(8, 6, [plotNetDir 'SingleNeuronActLevelFactorTimeCellType_' fileName])
    
    figure;
    hold on
    for nType = 1:3
        [fout, xout] = ksdensity(neuronOscTime(neuronType == nType), 0:0.02:1, 'Bandwidth', 0.02);
        stairs(xout, fout/sum(fout), 'linewid', 2)
        xlim([0 1])
    end
    box off
    xlabel('Oscillatory level')
    ylabel('Frac.')
    setPrint(8, 6, [plotNetDir 'SingleNeuronOscLevelFactorTimeCellType_' fileName])
    
    similarity_mat = nan(numNeuron, 2);
    for nNeuron = 1:numNeuron
        if ~isnan(neuronTime(nNeuron))
            nNeuronTime = neuronTime(nNeuron);
            if nNeuronTime < 10 || nNeuronTime > numTime - 75; continue; end
            dffs       = ca2spike(dff(nNeuron, :), new_activeNeuronMat(nNeuron, :), timePoints);
            % compute a minimal time window
            timeWindow = min(nNeuronTime, numTime - nNeuronTime - 40);
            x = dffs(1:timePoints(timeWindow-5));
            y = dffs((1:timePoints(timeWindow-5))+timePoints(nNeuronTime));
            z = dffs((1:timePoints(timeWindow-5))+timePoints(nNeuronTime + 40));
            [xcor, loc] = xcorr(x, y, 'coeff');
            similarity_mat(nNeuron, 1) = max(xcor(abs(loc)<1200));
            if isnan(similarity_mat(nNeuron, 1)); similarity_mat(nNeuron, 1) = 0; end
            [xcor, loc] = xcorr(z, y, 'coeff');
            similarity_mat(nNeuron, 2) = max(xcor(abs(loc)<1200));
        end
    end
    
    figure;
    subplot(1, 3, 1)
    hold on
    scatter(similarity_mat(:,1), similarity_mat(:,2), [], neuronType, 'filled')
    max_scale = ceil(nanmax(similarity_mat(:))*10)/10+0.01;
    xlim([0 max_scale])
    ylim([0 max_scale])
    refline(1)
    box off
    set(gca, 'TickDir', 'out');
    xlabel('Correlation before factored')
    ylabel('Correlation after factored')
    subplot(1, 3, 2)
    hold on
    for nType = 1:3
        if sum(~isnan(similarity_mat(neuronType == nType, 1))) ==0; continue; end
        [fout, xout] = ksdensity(similarity_mat(neuronType == nType, 1), 0:0.01:max_scale, 'Bandwidth', 0.01);
        stairs(xout, fout/sum(fout), 'linewid', 2)
        xlim([0 max_scale])
    end
    box off
    set(gca, 'TickDir', 'out');
    xlabel('Correlation before factored')
    ylabel('Frac.')
    subplot(1, 3, 3)
    hold on
    for nType = 1:3
        if sum(~isnan(similarity_mat(neuronType == nType, 1))) ==0; continue; end
        [fout, xout] = ksdensity(similarity_mat(neuronType == nType, 2), 0:0.01:max_scale, 'Bandwidth', 0.01);
        stairs(xout, fout/sum(fout), 'linewid', 2)
        xlim([0 max_scale])
    end
    box off
    set(gca, 'TickDir', 'out');
    xlabel('Correlation after factored')
    ylabel('Frac.')
    
    
    setPrint(8*3, 6, [plotNetDir 'SingleNeuronCorrLevelFactorTimeCellType_' fileName])
    
    save([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_1.mat'], 'neuronType', 'similarity_mat', '-append');
    
end

function dffs_marks = ca2spike(dffs, activeNeuronTime, timePoints)
    dffs_mark = ones(240, 1) * activeNeuronTime;
    dffs_mark = dffs_mark(:);
    [peaks, locs] = findpeaks(dffs(1:timePoints(numel(timePoints))), 'minPeakHeight', 0.0, 'minPeakDistance', 8);
    dffs_mark_locs = dffs_mark(locs);
    locs = locs(peaks > max([peaks(dffs_mark_locs==0), 0]) + 0.01);
    dffs_marks = false(size(dffs));
    dffs_marks(locs) = true;
    dffs_marks = double(dffs_marks);
end