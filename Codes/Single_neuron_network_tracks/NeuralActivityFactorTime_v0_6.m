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

function NeuralActivityFactorTime_v0_6(nFile)
    addpath('../Func');
    setDir;    
    fileName      = fileNames{nFile};   %#ok<*USENS>
    load([tempDatNetDir, fileName, '_spectrogram.mat'], 'peakPeriodMat');
    load([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_2.mat'], 'preLMat', 'preLMatIndex', 'preLMatTime') 
    numTime       = size(peakPeriodMat, 2);
    numFactor     = max(preLMatIndex);
    preLMatNeuron = sum(preLMat);
    timeWin       = 40;
    mColor        = cbrewer('qual', 'Dark2',  numFactor, 'cubic');
    nStep         = 20;
    
    peakPeriodMat(isnan(peakPeriodMat)) = 0;
    
    figure;
    hold on
    totNeuron     = 0;
    for nFactor   = 1:numFactor
        timeInd   = preLMatTime(:, preLMatNeuron<4 & preLMatIndex == nFactor);
        zeroTime  = min(timeInd);
        endTime   = min(max(timeInd), zeroTime+timeWin);
        LMat      = preLMat(:, preLMatNeuron<4 & preLMatIndex == nFactor & preLMatTime<=endTime);
        LMatInd   = sum(LMat, 2)>0;
        zeroTime  = min(timeInd);
        endTime   = min(max(timeInd), zeroTime+timeWin);
        minTime   = max(zeroTime - timeWin, 1);
        peakFrqValue  = bsxfun(@plus, peakPeriodMat(LMatInd, minTime:endTime), (1:sum(LMatInd))'*nStep + totNeuron*nStep);
        totNeuron = totNeuron + sum(LMatInd);
        hold on
        plot((minTime:endTime)-zeroTime, peakFrqValue, 'Color', mColor(nFactor,:), 'linewidth', 0.5);
        axis off
    end
    gridxy(0, [], 'color', 'w', 'linestyle', '--')
    xlim([-timeWin timeWin+1])
    ylim([-1 (totNeuron+1)*nStep+1])
%     setPrint(8, 6*length(LMatInd), [plotNetDir 'SingleNeuronDynamicsLocalCommunity_' fileName '_Factor_' num2str(nFactor, '%02d')], 'pdf')
end