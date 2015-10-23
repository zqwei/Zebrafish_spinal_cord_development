%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.  Covariance analysis: Factor analysis -- number of island using LONO 
% methods with selected neurons -- plots
% 
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function FACluster_v1_1(nFile)            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>    
    load([tempDatDir, fileName, '_numFactorNoLONOActiveNeuronsSimplified.mat']); 
    
    [numTime, numNeuron, ~] = size(EVLONO);
    
    figure;    
    [T, N] = meshgrid(1:numTime, 1:numNeuron);
    h = pcolor(T/60, N, squeeze(mean(EVLONO, 3)'));
    caxis([0 0.8])
    set(h, 'EdgeColor', 'none')
    axis xy
    xlabel('Time (hour)')
    ylabel('# latent dim.')
    xlim([0 numTime/60])
    set(gca, 'xtick', 1:6)
    colorbar
    box off
    setPrint(8, 6, [plotDir, 'FALONO_Time_Dim_', fileName], 'pdf')
    
    
    figure;    
    [T, N] = meshgrid(1:numTime, 1:numNeuron);
    [~, h] = contourf(T/60, N, squeeze(mean(EVLONO, 3)'));
    caxis([0 0.8])
    set(h, 'EdgeColor', 'none')
    axis xy
    xlabel('Time (hour)')
    ylabel('# latent dim.')
    xlim([0 numTime/60])
    set(gca, 'xtick', 1:6)
    colorbar
    box off
    setPrint(8, 6, [plotDir, 'FALONO_Time_Dim_Contour_', fileName], 'pdf')
    
end