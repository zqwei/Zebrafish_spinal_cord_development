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


function Neuron_selection_v3(nFile)
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

    activeNeuronMat  = false(size(dff, 1), length(timePoints));
    maxNeuronMat     = nan(size(dff, 1), length(timePoints));

%     lev = 5;
%     wname = 'sym8';
    
    for nNeuron    = 1:size(dff, 1)
        for nTime  = 1:length(timePoints)
            slicedDFF  = dff(nNeuron, timePoints(nTime)+1:timePoints(nTime)+timeStep);
%             powerNeuronMat(nNeuron, nTime)  = std(slicedDFF);
%             slicedDFF  = (slicedDFF - mean(slicedDFF))/std(slicedDFF);
            slicedDFF  = (slicedDFF - mean(slicedDFF))/std(slicedDFF);
%             [smoothedSlicedDFF, ~, ~, ~] = wden(slicedDFF, 'sqtwolog', 'h', 'mln', lev, wname);
            activeNeuronMat(nNeuron, nTime) = kstest2(-slicedDFF(slicedDFF<0), slicedDFF(slicedDFF>0), 'alpha', 0.05) && (skewness(slicedDFF)>0);
            maxNeuronMat(nNeuron, nTime)    = max(slicedDFF);
        end
    end
    
    refActiveNeuronMat = activeMatUseSGFit(rawf, 9, 511, timePoints, 0.05, timeStep);
    activeNeuronMat    = activeNeuronMat & refActiveNeuronMat;

%     for w                  = [21 41]
%         refActiveNeuronMat = activeMatUsePercentile(rawf, w, p, background, timePoints, 0.01);
%         activeNeuronMat    = activeNeuronMat & refActiveNeuronMat;
%     end
    
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
    
    save([tempDatDir, fileName, '.mat'], 'dff', 'tracks', 'leafOrder', 'slicedIndex', 'side', 'timePoints', 'sideSplitter', 'activeNeuronMat', 'timeStep');
    if exist('mnx', 'var')
        mnx       = mnx(slicedIndex); %#ok<NODEF>
        mnx       = mnx(leafOrder); %#ok<NASGU>
        save([tempDatDir, fileName, '.mat'], 'mnx', '-append')
    end
    
%     timeBin           = 11;
%     activeThres       = 5/timeBin;
%     
%     for nNeuron  = 1:size(activeNeuronMat, 1)
%         activeCurr = activeNeuronMat(nNeuron, :);
%         activeCurr = smooth(double(activeCurr), timeBin) > activeThres;
%         activeNeuronMat(nNeuron, timeBin:end) = activeCurr(timeBin:end);
%     end
%     
%     makeMovie(plotDir, fileName, timePoints, dff, activeNeuronMat, timeStep)
    
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

function makeMovie(plotDir, fileName, timePoints, dff, activeNeuronMat, timeStep)
%     if ispc
    video          = VideoWriter([plotDir '\movie_actMat_' fileName '.avi'], 'Uncompressed AVI');
%     elseif ismac
%         video          = VideoWriter([plotDir '\movie_actMat_' fileName '.mp4'], 'MPEG-4');
%     end
    video.FrameRate = 10;
    open(video);
    % frame size in pixels
    frameW = 1000;
    frameH = 1000;
    fig = figure('units', 'pixels', 'position', [0 0 , frameW, frameH]);
    set(0, 'defaultaxeslayer', 'top')
    linew = 1.25;
    for period = 1:numel(timePoints)
        timeRange = timePoints(period)+1:timePoints(period)+timeStep;
        activeTag = activeNeuronMat(:, period);
        clf reset
        hold on
        for i = 1:size(dff, 1)
            if ~activeTag(i)
                plot(linspace(0, 5, numel(timeRange)), zscore(dff(i, timeRange))+i*4, 'Color', 'b', 'linewidth', linew);
            else
                plot(linspace(0, 5, numel(timeRange)), zscore(dff(i, timeRange))+i*4, 'r', 'linewidth', linew);
            end
        end
        title(num2str(period))
        xlim([0, 5])
        ylim([0, size(dff, 1)*4+4]);
        set(gca, 'YTickLabel', '');
%         plot([0, 5], [sum(side==1) * 4, sum(side==1) * 4], 'k--');
        hold off
        frame = getframe(fig);
        writeVideo(video, frame);
    end
    close(video);
    close;
end