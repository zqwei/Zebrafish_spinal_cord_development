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


function Neuron_selection_v1_short_win(nFile)
    % load data
    addpath('../Func');
    setDir;

    fileDirName   = fileDirNames{nFile}; %#ok<*USENS>
    fileName      = fileNames{nFile};

    dirImageData  = [fileDirName '/'];
    load ([dirImageData 'profile.mat'])

    side = zeros(size(y));
    side(y<0) = 1;
    side(y>0) = 2;
    dff           = profile_all;
    rawf          = profile_all;
%     background    = 90;
%     baseline      = dff;
%     nCells        = size(dff, 1);
%     for nNeuron   = 1:nCells
%           baseline      = sgolayfilt(dff, 9, 511, [], 2);
%     end
%     dff           = bsxfun(@rdivide,(dff - baseline), (mean(baseline, 2)-background));
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
    
    timeStart     = 360;
    timeStep      = 720;
    numT          = size(dff, 2);
    timeEnd       = numT;
    
    slicedIndex   = ksTestAllTime;
    baseline      = baseline(slicedIndex, timeStart+1:timeEnd);
    rawf          = rawf(slicedIndex, timeStart+1:timeEnd);
    dff           = dff(slicedIndex, timeStart+1:timeEnd);
    tracks        = tracks(slicedIndex, timeStart+1:timeEnd, :);
    side          = side(slicedIndex);


    timePoints    = 0:240:timeEnd-timeStart-timeStep;
    slicedDFF       = dff;
    distNeurons     = pdist(slicedDFF, 'correlation');
    linkNeurons     = linkage(slicedDFF,'single','correlation');
    leafOrder       = optimalleaforder(linkNeurons, distNeurons);

    tSide         = side(leafOrder, :);
    sideIndex     = [find(tSide == 1); find(tSide == 2)];
    sideSplitter  = sum(tSide == 1)+0.5; %#ok<*NASGU>
    leafOrder     = leafOrder(sideIndex);

    baseline      = baseline(leafOrder, :);
    rawf          = rawf(leafOrder, :);
    dff           = dff(leafOrder, :); 
    tracks        = tracks(leafOrder, :, :);
    side          = side(leafOrder, :);

    activeNeuronMat   = false(size(dff, 1), length(timePoints));

    for nNeuron    = 1:size(dff, 1)
        for nTime  = 1:length(timePoints)
            slicedDFF           = dff(nNeuron, timePoints(nTime)+1:timePoints(nTime)+timeStep);
            slicedDFF           = (slicedDFF - mean(slicedDFF))/std(slicedDFF);
            activeNeuronMat(nNeuron, nTime) = kstest2(-slicedDFF(slicedDFF<0), slicedDFF(slicedDFF>0), 'alpha', 0.05) && (skewness(slicedDFF)>0);
        end
    end

%     makeMovie(plotDir, fileName, timePoints, dff, activeNeuronMat, timeStep)
    save([tempDatDir, fileName, '.mat'], 'dff', 'tracks', 'leafOrder', 'slicedIndex', 'side', 'timePoints', 'sideSplitter', 'activeNeuronMat', 'timeStep');
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
    video          = VideoWriter([plotDir '\movie_actMat_' fileName '.avi'], 'Uncompressed AVI');
    video.FrameRate = 10;
    open(video);
    frameW = 1000;
    frameH = 1000;
    fig = figure('units', 'pixels', 'position', [0 0 , frameW, frameH]);
    set(0, 'defaultaxeslayer', 'top')
    linew = 1.25;
    timeMax = timeStep/240;
    for period = 1:numel(timePoints)
        timeRange = timePoints(period)+1:timePoints(period)+timeStep;
        activeTag = activeNeuronMat(:, period);
        clf reset
        hold on
        for i = 1:size(dff, 1)
            if ~activeTag(i)
                plot(linspace(0, timeMax, numel(timeRange)), zscore(dff(i, timeRange))+i*4, 'Color', 'b', 'linewidth', linew);
            else
                plot(linspace(0, timeMax, numel(timeRange)), zscore(dff(i, timeRange))+i*4, 'r', 'linewidth', linew);
            end
        end
        title(num2str(period))
        xlim([0, timeMax])
        ylim([0, size(dff, 1)*4+4]);
        set(gca, 'YTickLabel', '');
        hold off
        frame = getframe(fig);
        writeVideo(video, frame);
    end
    close(video);
    close;
end