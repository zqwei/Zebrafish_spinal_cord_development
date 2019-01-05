% This code is to random crop series of trials from a long single trial
%
% version 1.0
%
%
% Input:
% SingleTrial :  Original single trial yDim x T (T -- time length)
% nLength     :  Time length for each cropped sections (nLength < T)
% nTrials     :  number of trials to crop
%
% Output:
% TrialSet    : yDim x nLength x nTrials

% Zqiang Wei
% Janelia Research Campus
% Email: weiz@janelia.hhmi.org

function TrialSet = randomCropSingleTrial(SingleTrial, nLength, nTrials)

    [yDim, T]      = size (SingleTrial);
    
    if nLength    >= T
        err('----- Trial length is smaller than that to crop -----');
    end
    
    TrialSet      = nan(yDim, nLength, nTrials);
    
    nStart        = floor(rand(nTrials,1) * (T - nLength -1)) + 1;
    
    for nTrial    = 1:nTrials
        TrialSet(:, :, nTrial)  = SingleTrial(:, nStart: nStart+nLength-1);
    end

end
