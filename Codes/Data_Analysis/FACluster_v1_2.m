%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.  Covariance analysis: Factor analysis -- loadin matrix evolution --
% reordering
% 
% size of factor against time
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function FACluster_v1_2(nFile)            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>    
    load([tempDatDir, 'EvoLoading_' fileName, '.mat'], 'networkMat')
    numTime           = length(networkMat);
    corrStd           = nan(numTime, 1);
    corrMean          = nan(numTime, 1);
    delayStd          = nan(numTime, 1);
    delayMean         = nan(numTime, 1);
    randomMean        = nan(numTime, 1);
    
    
    for nTime         = 1:numTime
        factorSet     = networkMat{nTime};
        dat           = [];
        for nFactor   = 2:length(factorSet)-1
            if length(factorSet{nFactor}.neuronIndex) > 1
                neuronIndex = factorSet{nFactor}.neuronIndex;
                neuronCCMat = factorSet{nFactor}.neuronCCMat;
                dat   = [dat; neuronCCMat(neuronIndex, :)]; %#ok<AGROW>
            end
            if size(dat, 1)>0
                datCorr           = dat(:,2);
                datCorr(isnan(datCorr)) = 0;
                corrStd(nTime)    = nanstd(datCorr);
                corrMean(nTime)   = nanmean(datCorr);
                delayStd(nTime)   = nanstd(abs(dat(:,1)));
                delayMean(nTime)  = nanmean(abs(dat(:,1)));
                randomMean(nTime) = mean(isnan(dat(:,1)));
            end
        end
    end
    
    figure;
    subplot(3, 1, 1)
    hold on
    plot((1:numTime)/60, corrMean, '-k', 'linewid', 2)
    plot((1:numTime)/60, corrMean+corrStd, '-k')
    plot((1:numTime)/60, corrMean-corrStd, '-k')
    hold off
    xlabel('Time (hour)')
    ylabel('Mean corr.')
    xlim([0 numTime/60])
    ylim([0 1])
    box off
    
    subplot(3, 1, 2)
    hold on
    plot((1:numTime)/60, delayMean, '-', 'linewid', 2, 'color', [     0    0.4470    0.7410])
    plot((1:numTime)/60, delayMean+delayStd, '-', 'color', [     0    0.4470    0.7410])
    plot((1:numTime)/60, delayMean-delayStd, '-', 'color', [     0    0.4470    0.7410])
    hold off
    xlabel('Time (hour)')
    ylabel('Mean delay (s)')
    xlim([0 numTime/60])
    ylim([0 25])
    box off
    
    subplot(3, 1, 3)
    plot((1:numTime)/60, randomMean, '-r', 'linewid', 2, 'color', [     0.8500    0.3250    0.0980])
    xlabel('Time (hour)')
    ylabel('Frac. of random delay')
    xlim([0 numTime/60])
    ylim([0 1])
    box off
    
    setPrint(8, 18, [plotDir, 'delay/FAIntraNeuron_', fileName], 'pdf');
end