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

    % Savitzky-Golay filtering
    sg_dff        = profile_all;
    numOrder      = 9;
    lenWindow     = 511;
    baseline      = sgolayfilt(sg_dff, numOrder, lenWindow, [], 2);
    sg_dff        = bsxfun(@rdivide,(sg_dff - baseline), mean(baseline, 2));
    nCells        = size(sg_dff, 1);

    % percentile window
    dff           = profile_all;
    w             = 20; % baselineWindowSize
    p             = 20; % baselinePrc
    baseline      = dff;
    for nNeuron   = 1:nCells
        baseline(nNeuron, :) = running_percentile(dff(nNeuron, :), w, p);
    end
    dff           = bsxfun(@rdivide,(dff - baseline), mean(baseline, 2));

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
    sg_dff        = sg_dff(slicedIndex, timeStep:timeEnd+timeStep);
    tracks        = tracks(slicedIndex, timeStep:timeEnd+timeStep, :);
    side          = side(slicedIndex); %#ok<NODEF>


    timePoints    = 0:240:timeEnd-timeStep;     %#ok<NASGU>


% %     logDetCorrAllDFF = zeros(length(timePoints),1);
% %     for nTime           = 1:length(timePoints)
% %         slicedDFF       = dff(:, timePoints(nTime)+1:timePoints(nTime)+1200);
% %         logDetCorrAllDFF(nTime) = log(det(corr(slicedDFF')));
% %     end
% %     figure;
% %     plot(timePoints/4/60, logDetCorrAllDFF, '-', 'linewid', 2)
% %     xlim([0 timePoints(end)/4/60]);
% %     xlabel('Time (min)')
% %     ylabel('log|Corr|')
% %     title('Colinearity index')
% %     box off
% %     setPrint(8, 6, [plotDir 'LogDetCorrIndex_' fileName], 'pdf')


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

% %     figure;
% %     for nPlot           = 1:10:length(timePoints)
% %         subplot(m, m, (nPlot-1)/10+1)
% %         hold on
% %         slicedDFF       = dff(:, timePoints(nPlot)+1:timePoints(nPlot)+1200);
% %         imagesc(corr(slicedDFF(leafOrder,:)'),[-1 1]);
% %         plot([1 numNeuron], [sideSplitter sideSplitter], '--w')
% %         plot([sideSplitter sideSplitter], [1 numNeuron], '--w')
% %         xlim([1 numNeuron])
% %         ylim([1 numNeuron])
% %         xlabel('Neuronal index')
% %         ylabel('Neuronal index')
% %         title([num2str(nPlot) 'min']);
% %         box off;
% %         hold off
% %     end
% %     setPrint(m*8, m*6, [plotDir, 'CorrMat_', fileName], 'pdf');
% %
% %
% %     figure;
% %     for nPlot           = 1:10:length(timePoints)
% %         subplot(m, m, (nPlot-1)/10+1)
% %         slicedDFF       = dff(:,timePoints(nPlot)+1:timePoints(nPlot)+1200);
% %         distNeurons     = pdist(slicedDFF, 'correlation');
% %         linkNeurons     = linkage(slicedDFF,'single','correlation');
% %         leafOrderLocal  = optimalleaforder(linkNeurons, distNeurons);
% %         imagesc(corr(slicedDFF(leafOrderLocal,:)'),[-1 1]);
% %         xlim([1 numNeuron])
% %         ylim([1 numNeuron])
% %         xlabel('Neuronal index')
% %         ylabel('Neuronal index')
% %         title([num2str(nPlot) 'min']);
% %         box off;
% %         axis xy
% %     end
% %     setPrint(m*8, m*6, [plotDir, 'CorrMatLocal_', fileName], 'pdf');


    baseline      = baseline(leafOrder, :); %#ok<NASGU>
    rawf          = rawf(leafOrder, :); %#ok<NASGU>
    dff           = dff(leafOrder, :); %#ok<NASGU>
    sg_dff        = sg_dff(leafOrder, :); %#ok<NASGU>
    tracks        = tracks(leafOrder, :, :); %#ok<NASGU>
    side          = side(leafOrder, :);  %#ok<NASGU>

    activeNeuronMat   = false(size(dff, 1), length(timePoints));
    activeNeuronSGMat = false(size(dff, 1), length(timePoints));

    for nNeuron    = 1:size(dff, 1)
        for nTime  = 1:length(timePoints)
            % Percentile activeNeuronMat
            slicedDFF           = dff(nNeuron, timePoints(nTime)+1:timePoints(nTime)+1200); %1200
            slicedDFF           = (slicedDFF - mean(slicedDFF))/std(slicedDFF);
            activeNeuronMat(nNeuron, nTime) = kstest2(-slicedDFF(slicedDFF<0), slicedDFF(slicedDFF>0), 'alpha', 0.01) && (skewness(slicedDFF)>0);
            % SG activeNeuronMat
            slicedDFF           = sg_dff(nNeuron, timePoints(nTime)+1:timePoints(nTime)+1200); %1200
            slicedDFF           = (slicedDFF - mean(slicedDFF))/std(slicedDFF);
            activeNeuronSGMat(nNeuron, nTime) = kstest2(-slicedDFF(slicedDFF<0), slicedDFF(slicedDFF>0), 'alpha', 0.05) && (skewness(slicedDFF)>0);
        end
    end

    activeNeuronMat = activeNeuronMat & activeNeuronSGMat;

    save([tempDatDir, fileName, '.mat'], 'dff', 'tracks', 'leafOrder', 'slicedIndex', 'side', 'timePoints', 'sideSplitter', 'activeNeuronMat');
    if exist('mnx', 'var')
        mnx       = mnx(slicedIndex); %#ok<NODEF>
        mnx       = mnx(leafOrder); %#ok<NASGU>
        save([tempDatDir, fileName, '.mat'], 'mnx', '-append')
    end

%     close all
%

end
