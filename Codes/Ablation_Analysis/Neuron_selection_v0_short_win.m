%
% Selecting neurons and storing their raw data, baseline, dff, tracks
%
% -------------------------------------------------------------------------
%
% version 0.0
% preprocessing of orginal data using window of r=10, prc=20%
% used for 3 min short win data (for ablation analysis)
% alpha=0.05 for before, alpha=0.01 for after 
%
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
%


function activeNeuronMat = Neuron_selection_v0_short_win(nFile)
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
    background    = 90;
    baseline      = dff;
    nCells        = size(dff, 1);
    w             = 21;
    p             = 20;
    
    for nNeuron   = 1:nCells
%           baseline      = sgolayfilt(dff, 9, 511, [], 2);
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
%     slicedDFF       = dff;
%     distNeurons     = pdist(slicedDFF, 'correlation');
%     linkNeurons     = linkage(slicedDFF,'single','correlation');
%     leafOrder       = optimalleaforder(linkNeurons, distNeurons);
    leafOrder     = 1:size(dff,1);

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

    alpha_value = 0.05;
    if mod(nFile, 2) == 0
        alpha_value = 0.01;
    end
    for nNeuron    = 1:size(dff, 1)
        for nTime  = 1:length(timePoints)
            slicedDFF           = dff(nNeuron, timePoints(nTime)+1:timePoints(nTime)+timeStep);
            slicedDFF           = (slicedDFF - mean(slicedDFF))/std(slicedDFF);
            activeNeuronMat(nNeuron, nTime) = kstest2(-slicedDFF(slicedDFF<0), slicedDFF(slicedDFF>0), 'alpha', 0.01) && (skewness(slicedDFF)>0);
        end
    end

    save([tempDatDir, fileName, '.mat'], 'dff', 'tracks', 'leafOrder', 'slicedIndex', 'side', 'timePoints', 'sideSplitter', 'activeNeuronMat', 'timeStep');
end

% function activeNeuronMat = activeMatUsePercentile(rawf, w, p, background, timePoints, alpha_value)
%     dff           = rawf;
%     baseline      = dff;
%     % percentile window
%     for nNeuron   = 1:size(dff, 1)
%         baseline(nNeuron, :) = running_percentile(dff(nNeuron, :), w, p);
%     end
%     dff           = bsxfun(@rdivide,(dff - baseline), (mean(baseline, 2)-background));
%     activeNeuronMat   = false(size(dff, 1), length(timePoints));
%     for nNeuron    = 1:size(dff, 1)
%         for nTime  = 1:length(timePoints)
%             slicedDFF           = dff(nNeuron, timePoints(nTime)+1:timePoints(nTime)+1200); %1200
%             slicedDFF           = (slicedDFF - mean(slicedDFF))/std(slicedDFF);
%             activeNeuronMat(nNeuron, nTime) = kstest2(-slicedDFF(slicedDFF<0), slicedDFF(slicedDFF>0), 'alpha', alpha_value) && (skewness(slicedDFF)>0);
%         end
%     end
% end
% 
% 
% function activeNeuronMat = activeMatUseSGFit(rawf, numOrder, lenWindow, timePoints, alpha_value)
%     dff           = rawf;
%     % Savitzky-Golay filtering
%     baseline      = sgolayfilt(dff, numOrder, lenWindow, [], 2);
%     dff           = bsxfun(@rdivide,(dff - baseline), mean(baseline, 2));
%     activeNeuronMat   = false(size(dff, 1), length(timePoints));
%     for nNeuron    = 1:size(dff, 1)
%         for nTime  = 1:length(timePoints)
%             slicedDFF           = dff(nNeuron, timePoints(nTime)+1:timePoints(nTime)+1200); %1200
%             slicedDFF           = (slicedDFF - mean(slicedDFF))/std(slicedDFF);
%             activeNeuronMat(nNeuron, nTime) = kstest2(-slicedDFF(slicedDFF<0), slicedDFF(slicedDFF>0), 'alpha', alpha_value) && (skewness(slicedDFF)>0);
%         end
%     end
% end

