%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  Covariance analysis: factor analysis
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.4 FA vs islands -- clustering -- Computing unrotated Loading matrix,
% Psi, correlation matrix etc. -- Raw data by FA communities -- KSTest for
% normality + spectrogram
%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function Data_Analysis_List2_7_7(nFile)

    if nargin<1
        nFile = 1;
    end
            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>
    
    load([tempDatDir, fileName, '.mat']);
    
    normalityMat      = zeros(size(dff, 1), length(timePoints)-1); %#ok<*NODEF>
    
    
    for nTime = 1:length(timePoints)-1
        slicedDFF     = dff(:,timePoints(nTime)+1:timePoints(nTime+1));
        slicedDFF     = bsxfun(@minus, slicedDFF, mean(slicedDFF,2));
        slicedDFF     = bsxfun(@rdivide, slicedDFF, std(slicedDFF,[],2))';
        for nNeuron   = 1:size(dff, 1);
            normalityMat(nNeuron, nTime) = kstest(slicedDFF(:, nNeuron), 'Alpha', 0.01);
        end
    end
    
    
%     figure;
%     imagesc(timePoints(1:end-1)/4/3600, 1:size(dff,1), 1 - normalityMat)
%     xlim([0 timePoints(end)'/4/3600])
%     ylim([1 size(dff,1)])
%     axis xy
%     colormap(gray)
%     colorbar
%     ylabel('Neuronal Index')
%     xlabel('Time (hour)')
%     title('Normality')  

    numNeuron     = size(dff, 1);
    m             = ceil(sqrt(numNeuron));
    
    window        = 1024;
    noverlap      = window/2;
    nfft          = window;
    fs            = 4;


    for nNeuron   = 1:numNeuron
        [~, f, t, pxx] = spectrogram(dff(nNeuron, :), window, noverlap, nfft, fs, 'power');
        subplot(m, m, nNeuron)
        f(f==0)   = 0.01;
        hold on;
        imagesc(t/3600, log2(f), log10(pxx));
        ylabel('log_2 Frequency (Hz)')
        xlabel('Time (hours)')
        axis xy
        plot(timePoints(1:end-1)/3600/4, normalityMat(nNeuron, :), '-k','linewid', 2)
        xlim([0 max(t)/3600])
    end    
    
    setFigureSize(8*m, 6*m)

end