%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.  Covariance analysis: Factor analysis -- number of island using non
% LONO methods with selected neurons
% 
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function FACluster_v0_2(nFile)    

    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>  
    load([tempDatDir, 'FALONO_', fileName, '.mat'], 'EVLONO');
    load([tempDatDir, fileName, '.mat'], 'dff', 'timePoints','activeNeuronMat');
    
    numTime           = size(EVLONO, 1);
    LONOM             = zeros(numTime, 1);
    
    EVLONOMat         = squeeze(mean(EVLONO, 3));
    
    for nTime         = 1:numTime    
        if sum(~isnan(EVLONOMat(nTime, :))) > 0
            maxEV        = nanmax(EVLONOMat(nTime, :));
            if maxEV>0
                LONOM(nTime) = find(EVLONOMat(nTime, :)>0.9*maxEV, 1, 'first');
            end
        end        
    end
    
    uncorrectedLONOM  = LONOM; %#ok<NASGU>
    
    maxNumFactor  = 12;
    numFactors    = repmat(struct('kgM', 0, 'paM', 0, 'SRMRM', 0, 'CFIM', 0), length(timePoints), 4); %#ok<NODEF>
    for nTime     = 1:length(timePoints)
        nNeuronIndex  = activeNeuronMat(:, nTime); %#ok<NODEF>
        if sum(nNeuronIndex)>1
            disp(nTime)
            slicedDFF = dff(nNeuronIndex,timePoints(nTime)+1:timePoints(nTime)+1200); 
            slicedDFF = bsxfun(@minus, slicedDFF, mean(slicedDFF,2));
            slicedDFF = bsxfun(@rdivide, slicedDFF, std(slicedDFF,[],2))';
            KG_M      = sum(eig(corr(slicedDFF))>1);
            currMaxFactor = min(KG_M, maxNumFactor);
            if currMaxFactor > 0
                numFactors(nTime, 1) = numFactorsWithNoCrossValidationSimplied(slicedDFF, currMaxFactor); 
            end
        end
    end 
    
    
    kgM   = [numFactors(:, 1).kgM]; 
    paM   = [numFactors(:, 1).paM]; 
    SRMRM = [numFactors(:, 1).SRMRM]; 
    CFIM  = [numFactors(:, 1).CFIM]; 
    
    timePoints        = (1:numTime)';
    
    figure;
    subplot(3, 2, 1)
    plot(timePoints/60, LONOM,'ko')
    ylabel('# factors')
    xlabel('Time (hour)')
    set(gca, 'xtick', 1:6)
    xlim([0 numTime/60])
    ylim([0 8])
    title('Leave one neuron out')
    box off
    
    subplot(3, 2, 3)
    plot(timePoints/60, kgM,'s','color',[0         0.4470    0.7410])
    ylabel('# factors')
    xlabel('Time (hour)')
    set(gca, 'xtick', 1:6)
    xlim([0 numTime/60])
    title('K-G')
    box off
    
    subplot(3, 2, 4)
    plot(timePoints/60, paM,'o','color',[0.8500    0.3250    0.0980])
    ylabel('# factors')
    xlabel('Time (hour)')
    set(gca, 'xtick', 1:6)
    xlim([0 numTime/60])
    title('Parallel analysis')
    box off
    
    subplot(3, 2, 5)
    plot(timePoints/60, SRMRM,'x','color',[0.9290    0.6940    0.1250])
    ylabel('# factors')
    xlabel('Time (hour)')
    set(gca, 'xtick', 1:6)
    xlim([0 numTime/60])
    title('SRMR')
    box off
    
    subplot(3, 2, 6)
    plot(timePoints/60, CFIM,'+','color',[0.4940    0.1840    0.5560])
    ylabel('# factors')
    xlabel('Time (hour)')
    set(gca, 'xtick', 1:6)
    xlim([0 numTime/60])
    title('CFI')
    box off

    setPrint(8*2, 6*3, [plotDir, 'numFactorLONOActiveNeurons_', fileName], 'pdf');
    close all
    
    save([tempDatDir, 'FALONO_', fileName, '.mat'], 'uncorrectedLONOM', 'numFactors', '-append');
        
end