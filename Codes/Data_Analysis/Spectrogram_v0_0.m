%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  Covariance analysis: spectrogram
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.4 FA vs islands -- clustering -- Computing unrotated Loading matrix,
% Psi, correlation matrix etc. -- Raw data by FA communities -- KSTest for
% distribution symmetry + spectrogram
%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function Spectrogram_v0_0(nFile)
            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>    
    load([tempDatDir, fileName, '.mat']);
    neuronIndex   = find(slicedIndex);
    neuronIndex   = neuronIndex(leafOrder);
        
    numNeuron     = size(dff, 1);
    m             = ceil(sqrt(numNeuron));
    window        = 1200;
    nw            = 7;
    fs            = 4;
    figure;
    
    timeWin       = 1200;
    lengthFreq    = timeWin/2;
    spectrogramMatAll = zeros(numNeuron, length(timePoints), lengthFreq);
    
    for nNeuron   = 1:numNeuron
        nTime    = 1;

        [pxx, f]       = pmtm(dff(nNeuron, timePoints(nTime)+1:timePoints(nTime)+timeWin), nw, window,fs);
        f              = f(2:end);
        spectrogramMat = zeros(length(timePoints), length(f));
        spectrogramMat(nTime, :) = pxx(2:end)./sum(pxx(2:end)); 

        for nTime     = 2:length(timePoints)
            pxx       = pmtm(dff(nNeuron, timePoints(nTime)+1:timePoints(nTime)+timeWin), nw, window,fs);
            spectrogramMat(nTime, :) = pxx(2:end)./sum(pxx(2:end));         
        end    
        
        spectrogramMatAll(nNeuron, :, :) = spectrogramMat;
        
        %%%%% Plot        
        subplot(m, m, nNeuron)
        hold on;
        imagesc(timePoints/60/60, f, log10(spectrogramMat'));
        ylabel('Frequency (s)')
        xlabel('Time (hour)')
        axis xy
        xlim([0 max(timePoints)/60/60])
        ylim([0 0.2])
        v = caxis;
        caxis([-2.5 v(2)]);
        colormap(jet)
        title(['Neuron #' num2str(neuronIndex(nNeuron))])
    end
    
    
    save([tempDatDir, fileName, '_spectrogram.mat'], 'spectrogramMatAll', 'f', '-append');
    setPrint(m*8, m*6, [plotDir, 'MTMSpectrogram_', fileName], 'pdf');

end


%     normalityMat      = zeros(size(dff, 1), length(timePoints)-1); %#ok<*NODEF>
%     
%     for nTime = 1:length(timePoints)-1
%         slicedDFF     = dff(:, timePoints(nTime)+1:timePoints(nTime)+1200);
%         slicedDFF     = bsxfun(@minus, slicedDFF, mean(slicedDFF,2));
%         slicedDFF     = bsxfun(@rdivide, slicedDFF, std(slicedDFF,[],2))';
%         for nNeuron   = 1:size(dff, 1);
%             nSlicedDFF = slicedDFF(:, nNeuron);
%             normalityMat(nNeuron, nTime) = kstest2(-nSlicedDFF(nSlicedDFF<0), nSlicedDFF(nSlicedDFF>0), 'Alpha', 0.01);
%         end
%     end
%     
