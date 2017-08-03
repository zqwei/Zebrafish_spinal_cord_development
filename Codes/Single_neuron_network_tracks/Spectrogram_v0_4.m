%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Overlay the wavelet analysis with active neuron index
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function Spectrogram_v0_4(nFile)
            
    addpath('../Func');
    setDir;    
    fileName      = fileNames{nFile}; 
    load([tempDatDir, fileName, '.mat'], 'timePoints');
    load([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_1.mat'], 'new_activeNeuronMat') 
    load([tempDatNetDir, fileName, '_spectrogram.mat'], 'spectrogramMatAll');
    
    timeWin       = 1200;
    numTime       = size(new_activeNeuronMat, 2);
    numNeuron     = size(new_activeNeuronMat, 1);
    m             = ceil(sqrt(numNeuron));
    
    specNeuronMat = nan(size(new_activeNeuronMat));
    
    
    for nNeuron   = 1:numNeuron
        spectrogramMat = mean(squeeze(spectrogramMatAll(nNeuron, :, :)));
        new_activeNeuronMat(nNeuron, :) = smooth(new_activeNeuronMat(nNeuron, :), 15) > 0;
        for nTime = 1:numTime
            specNeuronMat(nNeuron, nTime) = mean(spectrogramMat(timePoints(nTime)+(1:timeWin)));
        end
    end

% %     figure;
% %     for nNeuron   = 1:numNeuron
% %         %%%%% Plot        
% %         subplot(m, m, nNeuron)
% %         hold on;
% %         plot(specNeuronMat(nNeuron, :));
% %         plot(new_activeNeuronMat(nNeuron, :)*max(specNeuronMat(nNeuron, :)))
% %         gridxy([], prctile(specNeuronMat(nNeuron, new_activeNeuronMat(nNeuron, :)==0), 85))
% %         ylabel('Ave. Power')
% %         xlabel('Time (hour)')
% %     end
% %     
% %     setPrint(m*8, m*6, [plotNetDir, 'MTMSpectrogram_', fileName], 'pdf');
    
    oscNeuronMat = nan(size(new_activeNeuronMat));
    
    for nNeuron   = 1:numNeuron
        thres     = prctile(specNeuronMat(nNeuron, new_activeNeuronMat(nNeuron, :)==0), 85);
        if isempty(thres) || thres < 0 || isnan(thres)
            thres = 0;
        end
        oscNeuronMat(nNeuron, :) = specNeuronMat(nNeuron, :) > thres;
    end
    
    save([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_1.mat'], 'oscNeuronMat', '-append') 
    
end