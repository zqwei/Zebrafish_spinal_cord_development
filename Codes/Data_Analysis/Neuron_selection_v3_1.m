%
% Selecting neurons and storing their raw data, baseline, dff, tracks
%
% -------------------------------------------------------------------------
%
% version 0.0
% preprocessing of orginal data using sgolayfilt
%
%
% --- Deleting the time points has unreasonbale colinearity
% --- Deleting the units without oscillatory activity (KS test)
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
%


function Neuron_selection_v3_1(nFile, thresTwichCor)
    % load data
    addpath('../Func');
    setDir;

    fileDirName   = fileDirNames{nFile}; %#ok<USENS>
    fileName      = fileNames{nFile}; %#ok<USENS>

    dirImageData  = [fileDirName '/'];
    load ([dirImageData 'profile.mat'])

    % percentile window
    dff           = profile_all;
    rawf          = profile_all;
    w             = 61; % baselineWindowSize
    p             = 20; % baselinePrc
    background    = 90;
    baseline      = dff;
    nCells        = size(dff, 1);
    for nNeuron   = 1:nCells
        baseline(nNeuron, :) = running_percentile(dff(nNeuron, :), w, p);
    end
    dff           = bsxfun(@rdivide,(dff - baseline), (mean(baseline, 2)-background));

    tracks        = tracks_smoothed;

    ksTestAllTime    = false(nCells, 1);

    for nNeuron      = 1:nCells
        dat          = dff(nNeuron, :);
        if sum(isnan(dat)) == 0
            ksTestAllTime(nNeuron) = true;
        end
    end

    slicedIndex   = ksTestAllTime;


    timeStep      = 1200;
    numT          = size(dff, 2);
    timeEnd       = (floor(numT/timeStep) - 2)*timeStep;

    baseline      = baseline(slicedIndex, timeStep:timeEnd+timeStep);
    rawf          = rawf(slicedIndex, timeStep:timeEnd+timeStep);
    dff           = dff(slicedIndex, timeStep:timeEnd+timeStep);
    tracks        = tracks(slicedIndex, timeStep:timeEnd+timeStep, :);
    side          = side(slicedIndex); %#ok<NODEF>


    timePoints    = 0:240:timeEnd-timeStep;     %#ok<NASGU>

    slicedDFF       = dff(:, end - timeStep*10:end);
    distNeurons     = pdist(slicedDFF, 'correlation');
    linkNeurons     = linkage(slicedDFF,'single','correlation');
    leafOrder       = optimalleaforder(linkNeurons, distNeurons);

    tSide         = side(leafOrder, :);
    sideIndex     = [find(tSide == 1); find(tSide == 2)];
    sideSplitter  = sum(tSide == 1)+0.5; %#ok<NASGU>
    leafOrder     = leafOrder(sideIndex);

    numNeuron     = sum(slicedIndex);
    m             = ceil(sqrt(length(timePoints)/10+1));


    baseline      = baseline(leafOrder, :); %#ok<NASGU>
    rawf          = rawf(leafOrder, :); %#ok<NASGU>
    dff           = dff(leafOrder, :); %#ok<NASGU>
    tracks        = tracks(leafOrder, :, :); %#ok<NASGU>
    side          = side(leafOrder, :);  %#ok<NASGU>
    remvoeTwitchTime = false(size(rawf, 2), 1);

    % remove time points with twitching behavioral -- global negative
    % signal with strong positive correlation
    timeTwitch        = 60;
    thresTwitchNeuron = 10; % percentile of neurons
    if thresTwichCor > 0
        for nTime     = 1:length(dff)/timeTwitch
            slicedDFF = dff(:, (nTime-1)*timeTwitch+1:nTime*timeTwitch);
            negDFF    = slicedDFF < 0;
            [corVal, pVal] = corr(negDFF');
            corVal(pVal>0.05) = 0;
            if sum(sum(corVal > thresTwichCor)>thresTwitchNeuron) > thresTwitchNeuron
                dff(:, (nTime-1)*timeTwitch+1:nTime*timeTwitch) = nan;
            end
        end
        remvoeTwitchTime          = sum(isnan(dff))>0;
        dff(:, remvoeTwitchTime)  = [];
        tracks(:, remvoeTwitchTime, :) = [];
        timePoints    = timePoints(1:end-2);
    end

    activeNeuronMat  = false(size(dff, 1), length(timePoints));
    maxNeuronMat     = nan(size(dff, 1), length(timePoints));

    for nNeuron    = 1:size(dff, 1)
        for nTime  = 1:length(timePoints)
            slicedDFF  = dff(nNeuron, timePoints(nTime)+1:timePoints(nTime)+timeStep);
            slicedDFF  = (slicedDFF - mean(slicedDFF))/std(slicedDFF);
            activeNeuronMat(nNeuron, nTime) = kstest2(-slicedDFF(slicedDFF<0), slicedDFF(slicedDFF>0), 'alpha', 0.05) && (skewness(slicedDFF)>0);
            maxNeuronMat(nNeuron, nTime)    = max(slicedDFF);
        end
    end

    refActiveNeuronMat = activeMatUseSGFit(rawf(:, ~remvoeTwitchTime), 9, 511, timePoints, 0.05, timeStep);
    activeNeuronMat    = activeNeuronMat & refActiveNeuronMat;
    numActiveNeuron    = sum(activeNeuronMat, 1);

    for nTime      = 1:length(timePoints)
        slicedDFF  = dff(:, timePoints(nTime)+1:timePoints(nTime)+timeStep);
        [corrVal, p_val] = corr(slicedDFF');
        corrMat    = sum(corrVal > 0.3 & p_val < 0.05)>0;
        if sum(activeNeuronMat(:, nTime)) > 0
            actMax     = maxNeuronMat(activeNeuronMat(:, nTime), nTime);
            actThres   = min(actMax) + 0.9 * (max(actMax) - min(actMax));
            actThres   = max(actThres, prctile(maxNeuronMat(:, nTime), 90));
            actThres   = min(actThres, 5);
            activeNeuronMat(:, nTime) = activeNeuronMat(:, nTime) | maxNeuronMat(:, nTime) > actThres;
        else
            activeNeuronMat(:, nTime) = maxNeuronMat(:, nTime) > 5;
        end
        activeNeuronMat(:, nTime) = activeNeuronMat(:, nTime) | corrMat';
    end
    
    numActiveNewNeuron   = sum(activeNeuronMat, 1) - numActiveNeuron;
    save([tempDatDir, 'numActiveNewNeuron', fileName, '.mat'], 'numActiveNewNeuron');


end

function activeNeuronMat = activeMatUsePercentile(rawf, w, p, background, timePoints, alpha_value)
    dff           = rawf;
    baseline      = dff;
    % percentile window
    for nNeuron   = 1:size(dff, 1)
        baseline(nNeuron, :) = running_percentile(dff(nNeuron, :), w, p);
    end
    dff           = bsxfun(@rdivide,(dff - baseline), (mean(baseline, 2)-background));
    activeNeuronMat   = false(size(dff, 1), length(timePoints));
    for nNeuron    = 1:size(dff, 1)
        for nTime  = 1:length(timePoints)
            slicedDFF           = dff(nNeuron, timePoints(nTime)+1:timePoints(nTime)+timeStep); %timeStep
            slicedDFF           = (slicedDFF - mean(slicedDFF))/std(slicedDFF);
            activeNeuronMat(nNeuron, nTime) = kstest2(-slicedDFF(slicedDFF<0), slicedDFF(slicedDFF>0), 'alpha', alpha_value) && (skewness(slicedDFF)>0);
        end
    end
end


function activeNeuronMat = activeMatUseSGFit(rawf, numOrder, lenWindow, timePoints, alpha_value, timeStep)
    dff           = rawf;
    % Savitzky-Golay filtering
    baseline      = sgolayfilt(dff, numOrder, lenWindow, [], 2);
    dff           = bsxfun(@rdivide,(dff - baseline), mean(baseline, 2));
    activeNeuronMat   = false(size(dff, 1), length(timePoints));
    for nNeuron    = 1:size(dff, 1)
        for nTime  = 1:length(timePoints)
            slicedDFF           = dff(nNeuron, timePoints(nTime)+1:timePoints(nTime)+timeStep); %timeStep
            slicedDFF           = (slicedDFF - mean(slicedDFF))/std(slicedDFF);
            activeNeuronMat(nNeuron, nTime) = kstest2(-slicedDFF(slicedDFF<0), slicedDFF(slicedDFF>0), 'alpha', alpha_value) && (skewness(slicedDFF)>0);
        end
    end
end