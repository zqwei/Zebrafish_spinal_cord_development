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


function FACluster_v1_5           
    addpath('../Func');
    setDir;    
    
    networkSummary    = cell(length(fileNames), 1); %#ok<USENS>
    
    for nFile        = 1:length(fileNames)    
        fileName          = fileNames{nFile}; 
        load([tempDatDir, fileName, '.mat'], 'side', 'activeNeuronMat'); 
        load([tempDatDir, 'EvoLoading_' fileName, '.mat'], 'networkMat')
        numTime           = length(networkMat); %#ok<USENS>        
        activeNeuronMat   = sum(activeNeuronMat,1);
        dat               = {[], [], [], [], [], [], [], []};
    
        for nTime         = 1:numTime
            if activeNeuronMat(nTime) > 1
                factorSet     = networkMat{nTime};
                delayMat      = factorSet{end}.delayMat;
                corrMat       = factorSet{end}.corrMat;
                sizeMat       = ones(length(factorSet)-2);
                sideMat       = ones(length(factorSet)-2);
                for nFactor   = 2:length(factorSet)-1
                    sizeMat(nFactor-1) = length(factorSet{nFactor}.neuronIndex);
                    sideMat(nFactor-1) = mode(side(factorSet{nFactor}.neuronIndex));
                end
                for nFactor   = 1:length(corrMat)
                    for mFactor = nFactor+1:length(corrMat)
                        if sizeMat(nFactor)==1 && sizeMat(mFactor)==1 % neuron-neuron
                            if sideMat(nFactor) == sideMat(mFactor) 
                                dat{1} = [dat{1}; delayMat(nFactor, mFactor), corrMat(nFactor, mFactor)];
                            else
                                dat{2} = [dat{2}; delayMat(nFactor, mFactor), corrMat(nFactor, mFactor)];
                            end
                        elseif sizeMat(nFactor)>1 && sizeMat(mFactor)>1 % FA-FA
                            if sideMat(nFactor) == sideMat(mFactor) 
                                dat{5} = [dat{5}; delayMat(nFactor, mFactor), corrMat(nFactor, mFactor)];
                            else
                                dat{6} = [dat{6}; delayMat(nFactor, mFactor), corrMat(nFactor, mFactor)];
                            end
                        else % FA-neuron
                            if sideMat(nFactor) == sideMat(mFactor) 
                                dat{3} = [dat{3}; delayMat(nFactor, mFactor), corrMat(nFactor, mFactor)];
                            else
                                dat{4} = [dat{4}; delayMat(nFactor, mFactor), corrMat(nFactor, mFactor)];
                            end
                        end
                    end
                end
            end
        end
    end
    

    
    for nDat  = 1:length(saveFileNames)
        corrStd           = nan(numTime, 1);
        corrMean          = nan(numTime, 1);
        delayStd          = nan(numTime, 1);
        delayMean         = nan(numTime, 1);
        randomMean        = nan(numTime, 1);
        
        for nTime         = 1:numTime
            dat           = networkSummary{nTime}{nDat};
            if size(dat, 1)>0
                delayStd(nTime)   = nanstd(dat(:,1));
                delayMean(nTime)  = nanmean(dat(:,1));
                datCorr           = dat(:,2);
                datCorr(isnan(datCorr)) = 0;
                corrStd(nTime)    = nanstd(datCorr);
                corrMean(nTime)   = nanmean(datCorr);
                randomMean(nTime) = mean(isnan(dat(:,1)));
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

        setPrint(8, 18, [plotDir, 'delay/' saveFileNames{nDat}, fileName], 'pdf');
    end
end