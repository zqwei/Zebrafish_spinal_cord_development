%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spectrogram analysis of single neuron using wavelet 'amor'
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function Spectrogram_v0_3(nFile)
            
    addpath('../Func');
    setDir;    
    fileName      = fileNames{nFile}; 
    load([tempDatDir, fileName, '.mat']);
        
    numNeuron     = size(dff, 1);
    m             = ceil(sqrt(numNeuron));
    [b,a]         = butter(12,0.2,'low');
    sample_rate   = 4;
    
    for nNeuron   = 1:numNeuron
        dffs      = filtfilt(b, a, dff(nNeuron, :));
        if nNeuron == 1
            [cfs,f] = cwt(dffs,'amor', sample_rate);
            spectrogramMatAll = nan(size(dff, 1), length(f), size(dff, 2));
        else
            [cfs,~] = cwt(dffs,'amor', sample_rate);
        end
        spectrogramMatAll(nNeuron, :, :) = abs(cfs);
    end
    f_ind = f < 0.5 & f > 1/30; % long period as 30 ms
    spectrogramMatAll = spectrogramMatAll(:, f_ind, :);
    f     = f(f_ind);
    save([tempDatNetDir, fileName, '_spectrogram.mat'], 'spectrogramMatAll', 'f', '-v7.3');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% check with plots
%     load([tempDatNetDir, fileName, '_spectrogram.mat'], 'spectrogramMatAll', 'f');
%     for nNeuron   = 1:numNeuron
%         spectrogramMat = squeeze(spectrogramMatAll(nNeuron, :, :));
%         %%%%% Plot        
%         subplot(m, m, nNeuron)
%         hold on;
%         imagesc((1:size(spectrogramMatAll, 3))/4/60/60, f, spectrogramMat);
%         ylabel('Frequency (s)')
%         xlabel('Time (hour)')
%         axis xy
%         xlim([0 size(spectrogramMatAll, 3)/4/60/60])
%         ylim([0 0.5])
%         max_spec   = max(spectrogramMat(:));
%         caxis([max_spec*0.1 max_spec]);
%     end

    %%% check with mean power plots
    for nNeuron   = 1:numNeuron
        spectrogramMat = mean(squeeze(spectrogramMatAll(nNeuron, :, :)));
        %%%%% Plot        
        subplot(m, m, nNeuron)
        hold on;
        plot((1:size(spectrogramMatAll, 3))/4/60/60, spectrogramMat);
        gridxy([], 0.1*max(spectrogramMat), 'linestyle', '--')
        ylabel('Ave. Power')
        xlabel('Time (hour)')
        xlim([0 size(spectrogramMatAll, 3)/4/60/60])
    end
    
    setPrint(m*8, m*6, [plotNetDir, 'MTMSpectrogram_', fileName], 'pdf');

end