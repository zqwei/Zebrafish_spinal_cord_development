%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  Covariance analysis: spectrogram -- white noise test -- Ljung-Box Q
% test
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


function Spectrogram_v1(nFile)            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>   
    load([tempDatDir, fileName, '.mat'], 'sideSplitter')
    load([tempDatDir, fileName, '_spectrogram.mat'], 'spectrogramMatAll');
    
    timeWindow        = 1200;
    fs                = 4;
%     f                 = (1:timeWindow/2)/(timeWindow/2)*(fs/2);
    thresh            = 1/(timeWindow/2) * 3;
    
    
    [numNeuron, numTime, ~] = size(spectrogramMatAll); %#ok<NODEF>
    peakPeriodMat     = nan(numNeuron, numTime);
    
    for nNeuron       = 1:numNeuron
        for nTime     = 1:numTime
            sig       = squeeze(spectrogramMatAll(nNeuron, nTime, :));
            [peakLoc, peakMag] = peakfinder(sig, (max(sig)-min(sig))*0.95, thresh, 1, false, true); 
%             peakLoc   = peakfinder(sig, (max(sig)-min(sig))*0.95, [], 1, 0, true);            
            if ~isempty(peakLoc)
                peakFreq = peakLoc/(timeWindow/2)*(fs/2);
                peakMag  = peakMag(peakFreq>1/50 & peakFreq<1/5);
                peakFreq = peakFreq(peakFreq>1/50 & peakFreq<1/5);
                if ~isempty(peakFreq)
                    [~, maxFreq] = max(peakMag);
                    peakFreq  = peakFreq(maxFreq);
                    peakPeriodMat(nNeuron, nTime) = 1/peakFreq;
                end
            end
        end
    end
    
    figure;
    hold on;
    [T, N] = meshgrid(1:numTime, 1:numNeuron);
    h = pcolor(T/60, N, peakPeriodMat);
    colormap(cbrewer('div', 'RdYlGn',32))
    caxis([5 25])
    set(h, 'EdgeColor', 'none')
    set(gca, 'xtick', 1:6)
    plot([0 numTime/60], [sideSplitter sideSplitter], '--k')
    hold off
    xlim([0 numTime/60])
    ylim([0 numNeuron])
    axis xy
    xlabel('Time (hour)')
    ylabel('Neuron index')
    colorbar
    
    save([tempDatDir, fileName, '_spectrogram.mat'], 'peakPeriodMat', '-append');
    setPrint(8, 6, [plotDir, 'MTMSpectrogramPeak_', fileName], 'pdf');

end
