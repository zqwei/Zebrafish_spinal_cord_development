%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.  Covariance analysis: Factor analysis -- number of island using LONO 
% methods with selected neurons -- num dim from LONO
% 
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function FACluster_v1_2(nFile)            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>    
    load([tempDatDir, fileName, '_numFactorNoLONOActiveNeuronsSimplified.mat']); 
    
    numTime           = size(EVLONO, 1);
    LONOM             = zeros(numTime, 1);
    timePoints        = (1:numTime)';
    
    EVLONOMat         = squeeze(mean(EVLONO, 3));
    EVLONOMat(EVLONOMat==0) = nan;
    
    for nTime         = 1:numTime    
        if sum(~isnan(EVLONOMat(nTime, :))) > 0
            maxEV     = nanmax(EVLONOMat(nTime, :));
            LONOM(nTime) = find(EVLONOMat(nTime, :)>0.9*maxEV, 1, 'first');
        end        
    end
    
    figure;
    subplot(2, 3, 1)
    plot(timePoints/60, LONOM,'o')
    ylabel('# factors')
    xlabel('Time (hour)')
    set(gca, 'xtick', 1:6)
    xlim([0 numTime/60])
    title('Leave one neuron out')
    box off
    
    subplot(2, 3, 2)
    plot(timePoints/60, kgM,'o')
    ylabel('# factors')
    xlabel('Time (hour)')
    set(gca, 'xtick', 1:6)
    xlim([0 numTime/60])
    title('K-G')
    box off
    
    subplot(2, 3, 3)
    plot(timePoints/60, paM,'o')
    ylabel('# factors')
    xlabel('Time (hour)')
    set(gca, 'xtick', 1:6)
    xlim([0 numTime/60])
    title('Parallel analysis')
    box off
    
    subplot(2, 3, 4)
    plot(timePoints/60, SRMRM,'o')
    ylabel('# factors')
    xlabel('Time (hour)')
    set(gca, 'xtick', 1:6)
    xlim([0 numTime/60])
    title('SRMR')
    box off
    
    subplot(2, 3, 5)
    plot(timePoints/60, CFIM,'o')
    ylabel('# factors')
    xlabel('Time (hour)')
    set(gca, 'xtick', 1:6)
    xlim([0 numTime/60])
    title('CFI')
    box off

    setPrint(8*3, 6*2, [plotDir, 'numFactorActiveNeurons_', fileName], 'pdf');
    save([tempDatDir, fileName, '_numFactorNoLONOActiveNeuronsSimplified.mat'], 'LONOM', '-append'); 


end