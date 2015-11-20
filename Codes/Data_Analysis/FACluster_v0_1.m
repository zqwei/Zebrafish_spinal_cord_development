%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.  Covariance analysis: Factor analysis -- number of island using non
% LONO methods with all neurons
%
% -- selecting the active neurons...
% 
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function FACluster_v0_1(nFile)            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>    
    load([tempDatDir, 'FALONO_', fileName, '.mat'], 'EVLONO');
    
    [numTime, numNeuron, ~] = size(EVLONO);
    
    figure;    
    [T, N] = meshgrid(1:numTime, 1:numNeuron);
    h = pcolor(T/60, N, squeeze(mean(EVLONO, 3)'));
%     caxis([0 0.8])
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
%     caxis([0 0.8])
    set(h, 'EdgeColor', 'none')
    axis xy
    xlabel('Time (hour)')
    ylabel('# latent dim.')
    xlim([0 numTime/60])
    set(gca, 'xtick', 1:6)
    colorbar
    box off
    setPrint(8, 6, [plotDir, 'FALONOA_Time_Dim_Contour_', fileName], 'pdf')
    
    close all
    
end