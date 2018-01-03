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


function Neuron_selection_v1(nFile)
    % load data
    addpath('../Func');
    setDir;
    fileName      = fileNames{nFile}; %#ok<USENS>
    load([tempDatDir, fileName, '.mat'], 'dff', 'tracks', 'side', 'sideSplitter','timePoints', 'timeStep');


    numNeuron     = size(dff, 1); %#ok<NODEF>
    xTracks       = mean(mean(tracks, 3), 2);
    [~, leafOrder]= sortrows([xTracks, side], [2 -1]);
    m             = ceil(sqrt(length(timePoints)/10+1));

    figure;
    for nPlot           = 1:10:length(timePoints)
        subplot(m, m, (nPlot-1)/10+1)
        hold on
        slicedDFF       = dff(:, timePoints(nPlot)+1:timePoints(nPlot)+timeStep);
        % remove data with twitch times
        slicedDFF(:, sum(isnan(slicedDFF))>0)     = [];
        imagesc(corr(slicedDFF(leafOrder,:)'),[-1 1]);
        plot([1 numNeuron], [sideSplitter sideSplitter], '--w')
        plot([sideSplitter sideSplitter], [1 numNeuron], '--w')
        xlim([1 numNeuron])
        ylim([1 numNeuron])
        xlabel('Neuronal index')
        ylabel('Neuronal index')
        title([num2str(nPlot) 'min']);
        box off;
        hold off
    end
    setPrint(m*8, m*6, [plotDir, 'CorrMatLocation_', fileName], 'pdf');


end
