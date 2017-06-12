%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spectrogram analysis of single neuron using wavelet 'amor'
% This a plot from data in wavelet analysis
% data based on
%
% Spectrogram_v0_5
% 
% normalized power specturum by its local time peak
% maximum power specturm normalized by its max across time (in red)
% clustered time for each neuron, in white
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function Spectrogram_v0_5(nFile)
            
    addpath('../Func');
    setDir;    
    fileName      = fileNames{nFile}; 
    load([tempDatDir, fileName, '.mat']);
    
    load([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_1.mat'], 'new_activeNeuronMat') 
    load([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_2.mat'], 'preLMat', 'preLMatIndex', 'preLMatTime') 
        
    numNeuron     = size(dff, 1);
    m             = ceil(sqrt(numNeuron));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % check with plots
    load([tempDatNetDir, fileName, '_spectrogram.mat'], 'spectrogramMatAll', 'f');
    
    % for each neuron we define a matrix
    spec_mat_neuron = nan(length(f), length(timePoints));
    
    for nNeuron   = 1:numNeuron
        spectrogramMat = squeeze(spectrogramMatAll(nNeuron, :, :));
        
        for nTime = 1:length(timePoints)
            spec_mat_neuron(:, nTime) = mean(spectrogramMat(:, timePoints(nTime)+(1:1200)), 2);
        end
        
        ave_spec_mat                  = max(spec_mat_neuron, [], 1);
%         ave_spec_mat(ave_spec_mat<prctile(ave_spec_mat, 20)) = 1;
        ave_spec_mat(new_activeNeuronMat(nNeuron, :)==0) = 1;
        spec_mat_neuron_norm          = spec_mat_neuron./(ones(length(f),1) * ave_spec_mat);
        
        %%%% Plot        
%         subplot(m, m, nNeuron)
        ave_spec_mat                  = max(spec_mat_neuron, [], 1);
        ave_spec_mat                  = ave_spec_mat/max(ave_spec_mat);
        factorTime                    = preLMatTime(find(preLMat(nNeuron, :) == 1, 1, 'first'))/60;
        hold on;
        imagesc((timePoints/4/60+1)/60, f, spec_mat_neuron_norm);
        plot((timePoints/4/60+1)/60, ave_spec_mat * 0.5, '-r', 'linewid', 2)
        if ~isempty(factorTime)
            plot([factorTime factorTime], [0 0.5], '--w', 'linewid', 2)
        end
        ylabel('Frequency (s)')
        xlabel('Time (hour)')
        axis xy
        set(gca, 'TickDir', 'out')
        xlim([0 (timePoints(end)/4/60+2)/60])
        ylim([0 0.5])
        caxis([0 1])
        title(num2str(nNeuron))
    end
    
    setPrint(m*8, m*6, [plotNetDir, 'MTMSpectrogram_', fileName], 'pdf');

end