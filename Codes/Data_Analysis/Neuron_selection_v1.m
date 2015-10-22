%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0. Add side information as side index
% 
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function Neuron_selection_v1(nFile)            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>    
    load([tempDatDir, fileName, '.mat']); 
    load([tempDatDir, fileName, '_spectrogram.mat']);
    load([tempDatDir, fileName, '_numFactorNoLONOActiveNeuronsSimplified.mat']); 
    
    sideIndex         = [find(side == 1); find(side == 2)];
    baseline          = baseline(sideIndex, :);
    rawf              = rawf(sideIndex, :);
    dff               = dff(sideIndex, :);
    leafOrder         = leafOrder(sideIndex);
    tracks            = tracks(sideIndex, :);
    sideSplitter      = sum(side == 1)+0.5;
    side              = side(sideIndex, :);
    save([tempDatDir, fileName, '.mat'],...
        'dff', 'tracks', 'leafOrder', ...
        'baseline', 'rawf', ...
        'side', 'sideSplitter', '-append');
    
    figure;
    plot(side)
    
    spectrogramMatAll = spectrogramMatAll(sideIndex, :, :);
    peakPeriodMat     = peakPeriodMat(sideIndex, :);
    save([tempDatDir, fileName, '_spectrogram.mat'], 'spectrogramMatAll' , 'peakPeriodMat', '-append');
    
    sigNeuronsMat     = sigNeuronsMat(sideIndex, :);
    EVLONO(EVLONO == 0) = nan;
    
    save([tempDatDir, fileName, '_numFactorNoLONOActiveNeuronsSimplified.mat'],...
        'sigNeuronsMat', 'EVLONO', '-append');

end