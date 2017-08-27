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


function Neuron_selection_v0(nFile)
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

    activeNeuronMat   = false(size(dff, 1), length(timePoints));

    for nNeuron    = 1:size(dff, 1)
        for nTime  = 1:length(timePoints)
            slicedDFF           = dff(nNeuron, timePoints(nTime)+1:timePoints(nTime)+1200); %1200
            slicedDFF           = (slicedDFF - mean(slicedDFF))/std(slicedDFF);
            activeNeuronMat(nNeuron, nTime) = kstest2(-slicedDFF(slicedDFF<0), slicedDFF(slicedDFF>0), 'alpha', 0.01) && (skewness(slicedDFF)>0);
        end
    end

%     refActiveNeuronMat = activeMatUseSGFit(rawf, 9, 511, timePoints, 0.05);
%     activeNeuronMat    = activeNeuronMat & refActiveNeuronMat;
%     
%     for w                  = [21 41]
%         refActiveNeuronMat = activeMatUsePercentile(rawf, w, p, background, timePoints, 0.01);
%         activeNeuronMat    = activeNeuronMat & refActiveNeuronMat;
%     end

    save([tempDatDir, fileName, '.mat'], 'dff', 'tracks', 'leafOrder', 'slicedIndex', 'side', 'timePoints', 'sideSplitter', 'activeNeuronMat');
    if exist('mnx', 'var')
        mnx       = mnx(slicedIndex); %#ok<NODEF>
        mnx       = mnx(leafOrder); %#ok<NASGU>
        save([tempDatDir, fileName, '.mat'], 'mnx', '-append')
    end

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
            slicedDFF           = dff(nNeuron, timePoints(nTime)+1:timePoints(nTime)+1200); %1200
            slicedDFF           = (slicedDFF - mean(slicedDFF))/std(slicedDFF);
            activeNeuronMat(nNeuron, nTime) = kstest2(-slicedDFF(slicedDFF<0), slicedDFF(slicedDFF>0), 'alpha', alpha_value) && (skewness(slicedDFF)>0);
        end
    end
end


function activeNeuronMat = activeMatUseSGFit(rawf, numOrder, lenWindow, timePoints, alpha_value)
    dff           = rawf;
    % Savitzky-Golay filtering
    baseline      = sgolayfilt(dff, numOrder, lenWindow, [], 2);
    dff           = bsxfun(@rdivide,(dff - baseline), mean(baseline, 2));
    activeNeuronMat   = false(size(dff, 1), length(timePoints));
    for nNeuron    = 1:size(dff, 1)
        for nTime  = 1:length(timePoints)
            slicedDFF           = dff(nNeuron, timePoints(nTime)+1:timePoints(nTime)+1200); %1200
            slicedDFF           = (slicedDFF - mean(slicedDFF))/std(slicedDFF);
            activeNeuronMat(nNeuron, nTime) = kstest2(-slicedDFF(slicedDFF<0), slicedDFF(slicedDFF>0), 'alpha', alpha_value) && (skewness(slicedDFF)>0);
        end
    end
end

