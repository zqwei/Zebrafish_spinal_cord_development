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


function FACluster_v0_2_1(nFile)    

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
    
    LONOMLow          = zeros(numTime, 1);    
    EVLONOMat         = squeeze(mean(EVLONO, 3));
    for nTime         = 1:numTime    
        if sum(~isnan(EVLONOMat(nTime, :))) > 0
            maxEV        = nanmax(EVLONOMat(nTime, :));
            if maxEV>0
                LONOMLow(nTime) = find(EVLONOMat(nTime, :)>0.85*maxEV, 1, 'first');
            end
        end        
    end
    
    LONOMUp           = zeros(numTime, 1);    
    EVLONOMat         = squeeze(mean(EVLONO, 3));
    for nTime         = 1:numTime    
        if sum(~isnan(EVLONOMat(nTime, :))) > 0
            maxEV        = nanmax(EVLONOMat(nTime, :));
            if maxEV>0
                LONOMUp(nTime) = find(EVLONOMat(nTime, :)>0.95*maxEV, 1, 'first');
            end
        end        
    end
    
    timePoints        = (1:numTime)';
    figure;
    hold on    
    patch([timePoints'/60, flipud(timePoints/60)'], [LONOMLow', flipud(LONOMUp)'], 1,'facecolor',[0.5 0.5 0.5],'edgecolor','none');
%     alpha(h1, 0.5)
    plot(timePoints/60, LONOM,'ko')
    ylabel('# factors')
    xlabel('Time (hour)')
    set(gca, 'xtick', 1:6)
    xlim([0 numTime/60])
    ylim([0 8])
%     title('Leave one neuron out')
    box off


    setPrint(8, 6, [plotDir, 'numFactorLONOActiveNeuronsCI_', fileName]);
    
    close all
        
end