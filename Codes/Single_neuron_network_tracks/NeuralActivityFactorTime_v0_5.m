%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% draw neural activity for each small factor until the time small factor
% exploded or disappeared
% 
% computation from PreLMatTracker_v_0_1
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
    fileName      = fileNames{nFile};   %#ok<*USENS>
    load([tempDatNetDir, fileName, '_spectrogram.mat'], 'spectrogramMatAll', 'f');
    load([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_2.mat'], 'preLMat', 'preLMatIndex', 'preLMatTime') 
    numTime       = size(spectrogramMatAll, 2);
    numFactor     = max(preLMatIndex);
    preLMatNeuron = sum(preLMat);
    timeWin       = 40;
    mColor        = cbrewer('qual', 'Dark2',  numFactor, 'cubic');
    nStep         = 10;
            
    for nFactor   = 1:numFactor
        timeInd   = preLMatTime(:, preLMatNeuron<4 & preLMatIndex == nFactor);
        zeroTime  = min(timeInd);
        endTime   = min(max(timeInd), zeroTime+timeWin);
        LMat      = preLMat(:, preLMatNeuron<4 & preLMatIndex == nFactor & preLMatTime<=endTime);
        LMatInd   = find(sum(LMat, 2)>0);
        zeroTime  = min(timeInd);
        endTime   = min(max(timeInd), zeroTime+timeWin);
        minTime   = max(zeroTime - timeWin, 1);
        figure;
        max_spec  = max(max(max(spectrogramMatAll(LMatInd, minTime:endTime, :))));
        max_spec  = log10(max_spec);
        for mNeuron = 1:length(LMatInd)
            nNeuron = LMatInd(mNeuron);
            spectro = squeeze(spectrogramMatAll(nNeuron, minTime:endTime, :));
            subplot(length(LMatInd), 1, mNeuron)
            hold on;
            imagesc((minTime:endTime)-zeroTime, f, log10(spectro'));
            gridxy(0, [], 'color', 'w', 'linestyle', '--')
            ylabel('Frequency (s)')
            xlabel('Time (min)')
            axis xy
            xlim([-timeWin timeWin+1])
            ylim([0 0.2])
            caxis([-2.5 max_spec]);
            colormap(jet)
        end
        setPrint(8, 6*length(LMatInd), [plotNetDir 'SingleNeuronDynamicsLocalCommunity_' fileName '_Factor_' num2str(nFactor, '%02d')], 'pdf')
    end
end