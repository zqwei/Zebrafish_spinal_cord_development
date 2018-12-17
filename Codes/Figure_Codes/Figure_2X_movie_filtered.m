%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2X plot the evolution of factor center, together with atlas
%
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
addpath('../Func')
setDir;

nFile = 3;
tagNoLabel = 0;

fileName = fileNames{nFile};

load([tempDatDir, fileName, '.mat'], 'side', 'tracks', 'timePoints', 'activeNeuronMat', 'new_x', 'new_y', 'new_z', 'timeStep');
load([tempDatDir, 'LONOLoading_' fileName, '.mat'], 'CorrectedLMat')
    load([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_2.mat'], 'preLMat', 'preLMatIndex', 'preLMatTime') 


if tagNoLabel
    video          = VideoWriter([plotDir '\Movie_S7_smoothed_noLabel.avi'], 'Uncompressed AVI');
else
    video          = VideoWriter([plotDir '\Movie_S7_smoothed.avi'], 'Uncompressed AVI');
end

video.FrameRate = 10;
open(video);

if ~exist('new_x', 'var'); return; end

x                 = new_x;
y                 = new_y;
z                 = new_z;
numTime           = length(CorrectedLMat);
numNeuron         = length(side);



[~, neworder]     = sort(x);
neworder = [neworder(side(neworder)==1); neworder(side(neworder)==2)];

mColor = cbrewer('qual', 'Dark2',  8, 'cubic');
mColor = mColor([1 3 2 4:8], :);
mColor            = [mColor; cbrewer('qual', 'Set3',  15, 'cubic')];

linew = 1.25;
msScale = 3;
msOffset = 20;
radius = 0.2;


z = z/max(z) * 1.8;
y = y/2;

% frame size in pixels
frameW = 1200;
frameH = 650;
fig = figure('units', 'pixels', 'position', [0 0 , frameW, frameH]);
set(0, 'defaultaxeslayer', 'top')
whitebg('black');
set(gcf, 'Color', 'black');
set(0 ,'defaultfigurecolor', 'black')

subplot(1, 4, 1:3)
set(gca, 'ydir', 'reverse');
set(gca, 'xTick', 18:22);
xlim([17.5, 22]);
ylim([0 ceil(max(x))]);
set(gca, 'yTick', 0:ceil(max(x)));
if ~tagNoLabel
    xlabel('Age (hpf)');
    ylabel('AP location (segment)');
else
    axis off
end

factorIndexAll = preLMatIndex;
previousCenter = nan(size(factorIndexAll, 1), 4); %col#1: time point, XCenter, YCenter, flagUsed

for period = 1:numel(timePoints)
    LMat        = preLMat(:, preLMatTime==period);
    colorIndex = preLMatIndex(preLMatTime==period);
    activeTag = activeNeuronMat(:, period);

    
    timeRange = timePoints(period)+1:timePoints(period)+timeStep;
    subplot(1, 4, 4)
    cla reset
    set(gca, 'ydir', 'reverse');
    set(gca, 'xticklabel', '');
    set(gca, 'xtick', []);
    set(gca, 'yticklabel', '');
    set(gca, 'ytick', []);
    xlim([-1, 1])
    hold on
    
    if ~tagNoLabel
    box on
        for i = 1:ceil(max(x))-1
            plot([-2, 2], [i, i], '--w');
        end
        plot([0, 0], [0 ceil(max(x))], '--w');
    else
        axis off
    end
    plot(y(activeTag), x(activeTag), 'ow', 'MarkerFaceColor','w', 'linewidth', linew);
  
    for nFactor = 1:size(LMat, 2)
        % napolion plot
        subplot(1, 4, 1:3);
        hold on
        neuronFactor = LMat(:, nFactor)>0;
        dominateSide = ceil(mean(side(neuronFactor)) - 0.5);
        otherSide    = 3 - dominateSide;
        dominateSideNeuron = neuronFactor & side == dominateSide;
        xCenter = mean(x(dominateSideNeuron));
        yCenter = mean(y(dominateSideNeuron));
        scatter(timeRange(1)/(4*3600)+17.5, xCenter, sum(dominateSideNeuron)*10, mColor(colorIndex(nFactor), :), 'filled');
        
%         if ~isnan(previousCenter(colorIndex(nFactor), 1)) && previousCenter(colorIndex(nFactor), 4) == 0
%             plot([previousCenter(colorIndex(nFactor), 1), timeRange(1)/(4*3600)+17.5], [previousCenter(colorIndex(nFactor), 2), xCenter], 'Color',  mColor(colorIndex(nFactor), :));
%             previousCenter(colorIndex(nFactor), 4) = 1;
%         end
%         previousCenter(colorIndex(nFactor), 1) = timeRange(1)/(4*3600)+17.5;
%         previousCenter(colorIndex(nFactor), 2) = xCenter;
%         previousCenter(colorIndex(nFactor), 3) = yCenter;
%         previousCenter(colorIndex(nFactor), 4) = 0;
        % dorsal atlas
        subplot(1, 4, 4)       
        neuronFactor = LMat(:, nFactor)>0;
        dominateSide = ceil(mean(side(neuronFactor)) - 0.5);
        otherSideNeuron    = neuronFactor & side == otherSide;
        otherSide    = 3 - dominateSide;
        CHPoints = smoothedBoundary(x(dominateSideNeuron), y(dominateSideNeuron), radius);
        patch(CHPoints(:,2), CHPoints(:,1), mColor(colorIndex(nFactor), :), 'facealpha', 0.6, 'edgecolor', 'none');
        if sum(otherSideNeuron)>0
            if sum(otherSideNeuron) > 0
                CHPoints = smoothedBoundary(y(otherSideNeuron), x(otherSideNeuron), radius);
                patch(CHPoints(:,2), CHPoints(:,1), mColor(colorIndex(nFactor), :), 'facealpha', 0.6, 'edgecolor', 'none');
            end
        end
    end
    plot(y(~activeTag), x(~activeTag), 'ow', 'linewidth', linew, 'MarkerFaceColor', 'k');
    hold off
    frame = getframe(fig);
    writeVideo(video, frame);
end

close(video);
close;
