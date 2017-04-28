%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Factor movie
% 
% based on code -- 
% Yinan's version of FACluster_v0_5
% 
% -------------------------------------------------------------------------
% 
% Ziqiang Wei
% weiz@janelia.hhmi.org
%
%


function FactorMovie_v_0_1(nFile)            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile};   %#ok<*USENS>
    load([tempDatDir, fileName, '.mat'], 'dff', 'sideSplitter', 'side', 'tracks', 'timePoints', 'activeNeuronMat', 'new_x', 'new_y', 'new_z'); 
    load([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_1.mat'], 'CorrectedLMat', 'new_activeNeuronMat', 'lifeTimeTable') 
    
    if ~exist('new_x', 'var'); return; end
    
    x                 = new_x;
    y                 = new_y;
    z                 = new_z;
    numTime           = length(CorrectedLMat);
    numNeuron         = length(side);
    halfActTime       = nan(numNeuron, 1);
    timeBin           = 15;
    activeThres       = 0.65;
    lifeTimeThres     = 0.5;
    LMatThres         = 0.5;
    preLMat           = nan(numNeuron, 1);
    
    [~, neworder]     = sort(x);
    neworder = [neworder(side(neworder)==1); neworder(side(neworder)==2)];

    mColor = cbrewer('qual', 'Dark2',  8, 'cubic');
    mColor            = [mColor; cbrewer('qual', 'Set2',  128, 'cubic'); cbrewer('qual', 'Set2',  128, 'cubic'); cbrewer('qual', 'Set2',  128, 'cubic'); cbrewer('qual', 'Set2',  128, 'cubic')];
    
    
    if ispc
        video          = VideoWriter([plotNetDir '\movie_' fileName '.avi'], 'Uncompressed AVI');
    elseif ismac
        video          = VideoWriter([plotNetDir '\movie_' fileName '.mp4'], 'MPEG-4');
    end
    video.FrameRate = 10;
    open(video);


    % frame size in pixels
    frameW = 2000;
    frameH = 1000;
    fig = figure('units', 'pixels', 'position', [0 0 , frameW, frameH]);
    set(0, 'defaultaxeslayer', 'top')

    linew = 1.25;
    z = z/max(z) * 1.8;
    y = y/2;
    for period = 1:numel(timePoints)
        timeRange = timePoints(period)+1:timePoints(period)+1200;
        clf reset
        radius = 0.2;

        LMat      = CorrectedLMat{period};
        activeTag = new_activeNeuronMat(:, period);   %#ok<NODEF>
        lifeTime  = lifeTimeTable{period};
        
        % remove the factors with short life time
        LMat      = LMat(:, lifeTime >= lifeTimeThres);
        LMat      = LMat > LMatThres;
        LMat(:,sum(LMat, 1)==0) = [];
        
        % determine the factor index -- color index part
        if size(LMat,2) >= 1
            if sum(~isnan(preLMat(:))) == 0
                factorIndex  = 1:size(LMat, 2);
                preLMat      = LMat;
            else

%                 [~, maxFactorPerNeuronIndex] = max(LMat(sum(LMat, 2)>0, :), [], 2);
%                 sideRemoveList  = histc(maxFactorPerNeuronIndex, 1:size(LMat, 2)) <2; % remove the factor has no dominate factors
% 
%                 LMat(:, sideRemoveList) = [];
                
                sizeLMat        = sum(LMat, 1);
                [~, indexLMat]  = sort(sizeLMat, 'descend');
                LMat            = LMat(:, indexLMat);
                factorIndex     = zeros(size(LMat, 2), 1);
                similarityScore = zeros(size(LMat, 2), size(preLMat, 2));
                for nFactor     = 1:size(LMat, 2)
                    if sum(isnan(LMat(:)))>0; keyboard();end
                    similarityScore(nFactor, :) = sum(bsxfun(@and, LMat(:, nFactor), preLMat));
                end
                [~, maxIndex]   = max(similarityScore, [], 1);
                for nFactor     = 1:size(LMat, 2)
                    nFacotrNumPreFactor = sum(maxIndex == nFactor);
                    switch nFacotrNumPreFactor
                        case 0
                            preLMat = [preLMat, LMat(:, nFactor)];
                            factorIndex(nFactor) = size(preLMat, 2);
                        case 1
                            if similarityScore(nFactor, maxIndex == nFactor) == 0
                                preLMat = [preLMat, LMat(:, nFactor)];
                                factorIndex(nFactor) = size(preLMat, 2);
                            else
                                preLMatIndex         = find(maxIndex == nFactor);
                                factorIndex(nFactor) = preLMatIndex;
                                preLMat(:, preLMatIndex) = preLMat(:, preLMatIndex) | LMat(:, nFactor);
                            end
                        otherwise
                            preLMatIndex         = find(maxIndex == nFactor);
                            [~, nFactorMaxIndex] = max(similarityScore(nFactor, preLMatIndex));
                            factorIndex(nFactor) = preLMatIndex(nFactorMaxIndex);
                            preLMat(:, preLMatIndex(nFactorMaxIndex)) = preLMat(:, preLMatIndex(nFactorMaxIndex)) | LMat(:, nFactor);
                    end
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plots
        
        % plot calcium traces
        subplot(3, 2, [1 3 5]);
        hold on
        for i = 1:size(LMat, 1)
            if ~activeTag(i)
                plot(linspace(0, 5, numel(timeRange)), zscore(dff(i, timeRange))+find(neworder==i)*4, 'Color', [.8, .8, .8], 'linewidth', linew);
            end
        end
        for i = 1:size(LMat, 1)
            if activeTag(i)
                if ~any(isnan(LMat(i, :))) && sum(LMat(i, :))>0 && size(LMat,2) >= 1
                    [~, nFactor] = max(LMat(i, :));
                    plot(linspace(0, 5, numel(timeRange)), zscore(dff(i, timeRange))+find(neworder==i)*4, 'Color', mColor(factorIndex(nFactor), :), 'linewidth', linew);
                else
                    plot(linspace(0, 5, numel(timeRange)), zscore(dff(i, timeRange))+find(neworder==i)*4, 'k', 'linewidth', linew);
                end
            end
        end
        xlim([0, 5])
        ylim([0, size(LMat, 1)*4+4]);
        set(gca, 'YTickLabel', '');
        %     text(zeros(numel(neworder), 1), (1:numel(neworder))*4, num2str(neworder));
        plot([0, 5], [sum(side==1) * 4, sum(side==1) * 4], 'k--');
        hold off

        % plot spatial organization
        % dorsal view
        subplot(3, 2, 2)
        hold on
        plot(x(activeTag), y(activeTag), 'ok', 'MarkerFaceColor','k', 'linewidth', linew);
        if size(LMat,2) >= 1
            for nFactor = 1:size(LMat, 2)
                neuronFactor = LMat(:, nFactor)>0;
                dominateSide = ceil(mean(side(neuronFactor)) - 0.5);
                otherSide    = 3 - dominateSide;
                dominateSideNeuron = neuronFactor & side == dominateSide;
                otherSideNeuron    = neuronFactor & side == otherSide;
                CHPoints = smoothedBoundary(x(dominateSideNeuron), y(dominateSideNeuron), radius);
                patch(CHPoints(:,1), CHPoints(:,2), mColor(factorIndex(nFactor), :), 'facealpha', 0.6, 'edgecolor', 'none');
                if sum(otherSideNeuron)>0
                    if sum(otherSideNeuron) > 0
                        CHPoints = smoothedBoundary(x(otherSideNeuron), y(otherSideNeuron), radius);
                        patch(CHPoints(:,1), CHPoints(:,2), mColor(factorIndex(nFactor), :), 'facealpha', 0.6, 'edgecolor', 'none');
                    end
                end
                
%                 if lifeTime(nFactor) < lifeTimeThres
%                     plot(x(neuronFactor), y(neuronFactor), 'or', 'MarkerFaceColor','r', 'linewidth', linew);
%                 end
                
            end
        end
        text(floor(max(x))-1, 0, num2str(period),'color','r','fontsize', 24)
        plot(x(~activeTag), y(~activeTag), 'ok', 'linewidth', linew, 'MarkerFaceColor', 'w');
        for i = 1:ceil(max(x))
            plot([i, i], [-2, 2], '--k');
        end
        xlim([0 ceil(max(x))+1]);
        ylim([-1 1]);
        hold off

        % side view - left
        subplot(3, 2, 6)
        hold on
        plot(x(activeTag & side==1), z(activeTag & side==1), 'ok', 'MarkerFaceColor','k', 'linewidth', linew);
        if size(LMat,2) >= 1
            for nFactor = 1:size(LMat, 2)
                neuronFactor = LMat(:, nFactor)>0;
                if sum(neuronFactor & side==1)>0
                    CHPoints = smoothedBoundary(x(neuronFactor & side==1), z(neuronFactor & side==1), radius);
                    patch(CHPoints(:,1), CHPoints(:,2), mColor(factorIndex(nFactor), :), 'facealpha', 0.6, 'edgecolor', 'none');
                end
            end
        end
        plot(x(~activeTag & side==1), z(~activeTag & side==1), 'ok', 'linewidth', linew, 'MarkerFaceColor', 'w');
        for i = 1:ceil(max(x))
            plot([i, i], [0, 2], '--k');
        end
        xlim([0 ceil(max(x))+1]);
        ylim([0 2]);
        hold off

        % side view - right
        subplot(3, 2, 4)
        hold on
        plot(x(activeTag & side==2), z(activeTag & side==2), 'ok', 'MarkerFaceColor','k', 'linewidth', linew);
        if size(LMat,2) >= 1
            for nFactor = 1:size(LMat, 2)
                neuronFactor = LMat(:, nFactor)>0;
                if sum(neuronFactor & side==2)>0
                    CHPoints = smoothedBoundary(x(neuronFactor & side==2), z(neuronFactor & side==2), radius);
                    patch(CHPoints(:,1), CHPoints(:,2), mColor(factorIndex(nFactor), :), 'facealpha', 0.6, 'edgecolor', 'none');
                end
            end
        end
        plot(x(~activeTag & side==2), z(~activeTag & side==2), 'ok', 'linewidth', linew, 'MarkerFaceColor', 'w');
        for i = 1:ceil(max(x))
            plot([i, i], [0, 2], '--k');
        end
        xlim([0 ceil(max(x))+1]);
        ylim([0 2]);
        hold off
        frame = getframe(fig);
        writeVideo(video, frame);
    end
    close(video);
    close;
    
end