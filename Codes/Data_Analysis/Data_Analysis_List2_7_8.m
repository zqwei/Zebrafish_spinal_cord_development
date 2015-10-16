%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  Covariance analysis: factor analysis
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.4 FA vs islands -- clustering -- Computing unrotated Loading matrix,
% Psi, correlation matrix etc. -- Raw data by FA communities -- FFT --
% Welch's methods
%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function Data_Analysis_List2_7_8(nFile)

    if nargin<1
        nFile = 1;
    end
    
    addpath('../Func');
    setDir;    
    fileName      = fileNames{nFile}; %#ok<USENS>
    
    load([tempDatDir, fileName, '.mat'], 'dff','timePoints');
    
    numNeuron     = size(dff, 1); %#ok<NODEF>
        
    PSDPeakTime   = nan(size(dff, 1), length(timePoints)-1);
    hannWindow    = 1024;
    noverlap      = hannWindow/2;
    nfft          = hannWindow;
    fs            = 4;
    vecHannWin    = hann(hannWindow);
    
    numOrder           = 19;
    lenWindow          = 511;
    
    for nTime = 1:length(timePoints)-1
        slicedDFF     = dff(:,timePoints(nTime)+1:timePoints(nTime+1));
        slicedDFF     = bsxfun(@minus, slicedDFF, mean(slicedDFF,2));
        slicedDFF     = bsxfun(@rdivide, slicedDFF, std(slicedDFF,[],2))';
        for nNeuron   = 1:numNeuron
            if kstest(slicedDFF(:, nNeuron))
                nSlicedDFF       = slicedDFF(:, nNeuron);
                stdNSlicedDFF    = (nSlicedDFF - mean(nSlicedDFF))/std(nSlicedDFF);
                filtedNSlicedDFF = stdNSlicedDFF - sgolayfilt(stdNSlicedDFF, numOrder, lenWindow);
                [fftDff, f] = periodogram(filtedNSlicedDFF, [], nfft, fs, 'power');
%                 [fftDff, f] = pwelch(filtedNSlicedDFF, vecHannWin, noverlap, nfft, fs);
                [~, maxFIndex] = max(fftDff);
                PSDPeakTime(nNeuron, nTime) = f(maxFIndex);
            end
        end
    end
    
    figure;
    h = imagesc(timePoints(1:end-1)/4/3600, 1:size(dff,1), 1./PSDPeakTime); %#ok<COLND>
    xlim([0 timePoints(end)'/4/3600]) %#ok<COLND>
    ylim([1 size(dff,1)])
    ylabel('Neuronal Index')
    xlabel('Time (hour)')
    title('Power Spectrum Peak Cycle Width')
    colormap(linspecer(nfft/fs)), colorbar
    set(h,'AlphaData',~isnan(PSDPeakTime) & ~isinf(1./PSDPeakTime))
    
    save([tempDatDir, fileName, '_PSDPeakTime.mat'], 'PSDPeakTime'); 
    
%% Computation through other methods (sample frequency nfft is also tuned)
%
%     FFTPeakTime   = zeros(size(dff, 1), length(timePoints)-1);
%     FFTMTMPeakTime = zeros(size(dff, 1), length(timePoints)-1);
%     PSDMinHoldPeakTime = zeros(size(dff, 1), length(timePoints)-1);
% 
%     Fs            = 4; % sample rate
%     TW            = 3; % taper time-bandwidth
%     params.tapers = [TW 2*TW-1];
%     params.pad    = -1;
%     params.Fs     = Fs;
% %     pErr          = 0.05;
% %     params.err    = [1 pErr];
% 
%     
%     
%     for nTime = 1:length(timePoints)-1
%         slicedDFF     = dff(:,timePoints(nTime)+1:timePoints(nTime+1));
%         slicedDFF     = bsxfun(@minus, slicedDFF, mean(slicedDFF,2));
%         slicedDFF     = bsxfun(@rdivide, slicedDFF, std(slicedDFF,[],2))';
%         
%         %%% FFT results
%         lengthDff     = size(slicedDFF, 1);
%         f             = Fs*(0:(lengthDff/2))/lengthDff;
%         fftDff        = abs(fft(slicedDFF, [], 1)/lengthDff);
%         % only consider those in between 5s to 100s cycles
%         fftDff        = fftDff(f>=0.01 & f<=0.2, :);
%         f             = f(f>=0.01 & f<=0.2);
%         [~, fftPeakIndex] = max(fftDff, [], 1);
%         FFTPeakTime(:, nTime) = f(fftPeakIndex);
%         
%         %%% FFT-MTM results
%         [S ,f]        = mtspectrumc(slicedDFF,params);
%         % only consider those in between 5s to 200s cycles
%         % also simplify the computation for fftDff = sqrt(S/(lengthDff/TW/Fs))
%         fftDff        = S(f>=0.005 & f<=0.2, :);
%         f             = f(f>=0.005 & f<=0.2);
%         [~, fftPeakIndex] = max(fftDff, [], 1);
%         FFTMTMPeakTime(:, nTime) = f(fftPeakIndex);
%         
%         [fftDff, f] = pwelch(slicedDFF, hann(800), 400, 800, 4);
%         fftDff      = bsxfun(@rdivide, fftDff, max(fftDff,[],1));
%         for nNeuron = 1:size(dff, 1) 
%             [fftPeaks, fftPeakIndex] = findpeaks(fftDff(:, nNeuron), f,'SortStr','descend','NPeaks',2);
%             fftPeaks  = fftPeaks(fftPeakIndex > 0.025 & fftPeakIndex < 1);
%             fftPeakIndex  = fftPeakIndex(fftPeakIndex > 0.025 & fftPeakIndex < 1);
%             if isempty(fftPeaks) || max(fftPeaks) < 0.1
%                 PSDMinHoldPeakTime(nNeuron, nTime) = nan;
%             else
%                 PSDMinHoldPeakTime(nNeuron, nTime) = fftPeakIndex(fftPeaks == max(fftPeaks));
%             end
%         end
%         
%     end
    
    
    
    
end
