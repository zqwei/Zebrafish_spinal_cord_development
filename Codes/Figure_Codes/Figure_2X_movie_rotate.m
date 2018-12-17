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


fileName = fileNames{nFile};

load([tempDatDir, fileName, '.mat'], 'side', 'tracks', 'timePoints', 'activeNeuronMat', 'new_x', 'new_y', 'new_z', 'timeStep');
load([tempDatDir, 'LONOLoading_' fileName, '.mat'], 'CorrectedLMat')

[~, neworder]     = sort(x);
neworder = [neworder(side(neworder)==1); neworder(side(neworder)==2)];

mColor = cbrewer('qual', 'Dark2',  8, 'cubic');
mColor            = [mColor; cbrewer('qual', 'Set2',  128, 'cubic')];
preLMat           = nan(numNeuron, 1);

video          = VideoWriter([plotDir '\Movie_S7_rotation.avi'], 'Uncompressed AVI');

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
mColor            = [mColor; cbrewer('qual', 'Set2',  128, 'cubic')];
preLMat           = nan(numNeuron, 1);

linew = 1.25;
msScale = 3;
msOffset = 20;

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
xlabel('Age (hpf)');
ylabel('AP location (segment)');

elSeries = linspace(-90, 90, (numel(timePoints)));

for period = 1:numel(timePoints)
    timeRange = timePoints(period)+1:timePoints(period)+timeStep;
    radius = 0.2;
    
    LMat          = CorrectedLMat{period};
    LMat(isnan(LMat)) = 0;
    LMat(:, sum(LMat, 1)==0) = [];
    activeTag = activeNeuronMat(:, period);
    
    % code to drop overlapped factors
    LMatNeuron    = LMat>0;
    numFactor     = size(LMatNeuron, 2);
    for nFactor   = 1:size(LMatNeuron, 2)
        if sum(LMatNeuron(:,nFactor)) >0 % skip if nFactor is dropped
            for mFactor = nFactor+1:size(LMatNeuron, 2)
                if sum(LMatNeuron(:,mFactor)) >0 % skip if mFactor is dropped
                    if all(ismember(find(LMatNeuron(:,nFactor))', find(LMatNeuron(:,mFactor))')) % if nFactor inside any other factors, drop nFactor
                        LMatNeuron(:,nFactor) = 0;
                        continue;
                    elseif all(ismember(find(LMatNeuron(:,mFactor))', find(LMatNeuron(:,nFactor))')) % if mFactor inside any other factors, drop nFactor
                        LMatNeuron(:,mFactor) = 0;
                        continue;
                    end
                end
            end
        end
    end
    LMat(LMatNeuron == 0)    = 0;
    LMat(:, sum(LMat, 1)==0) = []; % drop factors with zero weight
    % end of drop overlapped factor code
    
    % determine the factor index
    % break double sided factors
    LMatLeft                 = zeros(size(LMat));
    LMatRight                = zeros(size(LMat));
    LMatLeft(side == 1, :)   = LMat(side == 1, :);
    LMatRight(side == 2, :)  = LMat(side == 2, :);
    LMat                     = [LMatLeft, LMatRight];
    
    % remove the single-unit factor in movie
    LMat(:, sum(LMat>0, 1)<=1) = []; % drop factors with zero weight
    
    if size(LMat,2) >= 1
        if sum(~isnan(preLMat(:))) == 0
            factorIndex  = 1:size(LMat, 2);
            preLMat      = LMat;
        else
            [~, maxFactorPerNeuronIndex] = max(LMat(sum(LMat, 2)>0, :), [], 2);
            sideRemoveList  = histc(maxFactorPerNeuronIndex, 1:size(LMat, 2)) <2; % remove the factor has no dominate factors
            
            LMat(:, sideRemoveList) = [];
            LMat            = LMat > 0;
            sizeLMat        = sum(LMat, 1);
            [~, indexLMat]  = sort(sizeLMat, 'descend');
            LMat            = LMat(:, indexLMat);
            factorIndex     = zeros(size(LMat, 2), 1);
            
            % compute similarity matrix
            similarityScore = zeros(size(LMat, 2), size(preLMat, 2));
            for nFactor     = 1:size(LMat, 2)
                if sum(isnan(LMat(:)))>0; keyboard();end
                similarityScore(nFactor, :) = sum(bsxfun(@and, LMat(:, nFactor), preLMat));
            end
            
            % check if any prefactor has no connection with new factors
            % decide which factor is not included in preLMatIndex
            % maxIndex is the factor with the maximum coverage with the prefactors
            [~, maxIndex]   = max(similarityScore, [], 1);
            
            % check if any prefactor is merged (factor has maximum coverages with more than one prefactors, pick the larger one as its index)
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
        subplot(1, 4, 4)
        cla reset
        set(gca, 'ydir', 'reverse');
        set(gca, 'xticklabel', '');
        set(gca, 'xtick', []);
        set(gca, 'yticklabel', '');
        set(gca, 'ytick', []);
        xlim([-1, 1])
        box on
        hold on
        for i = 1:ceil(max(x))-1
            plot([-2, 2], [i, i], '--w');
        end
        plot([0, 0], [0 ceil(max(x))], '--w');
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
            %             scatter(timeRange(1)/(4*3600)+17.5, xCenter, sum(dominateSideNeuron)*10, mColor(factorIndex(nFactor), :), 'filled');
%             scatter(timeRange(1)/(4*3600)+17.5, xCenter, sum(dominateSideNeuron)*msScale + msOffset, mColor(factorIndex(nFactor), :), 'filled');
            scatter3(timeRange(1)/(4*3600)+17.5, xCenter, yCenter, sum(dominateSideNeuron)*msScale + msOffset, mColor(factorIndex(nFactor), :), 'filled');
            % dorsal atlas
            subplot(1, 4, 4)         
            

            neuronFactor = LMat(:, nFactor)>0;
            dominateSide = ceil(mean(side(neuronFactor)) - 0.5);
            otherSideNeuron    = neuronFactor & side == otherSide;
            otherSide    = 3 - dominateSide;
            CHPoints = smoothedBoundary(x(dominateSideNeuron), y(dominateSideNeuron), radius);
            patch(CHPoints(:,2), CHPoints(:,1), mColor(factorIndex(nFactor), :), 'facealpha', 0.6, 'edgecolor', 'none');
            if sum(otherSideNeuron)>0
                if sum(otherSideNeuron) > 0
                    CHPoints = smoothedBoundary(y(otherSideNeuron), x(otherSideNeuron), radius);
                    patch(CHPoints(:,2), CHPoints(:,1), mColor(factorIndex(nFactor), :), 'facealpha', 0.6, 'edgecolor', 'none');
                end
            end           
        end
        plot(y(~activeTag), x(~activeTag), 'ow', 'linewidth', linew, 'MarkerFaceColor', 'k');
        hold off
    end
    subplot(1, 4, 1:3) 
    view([0, elSeries(period)])
    frame = getframe(fig);
    writeVideo(video, frame);
end

close(video);
close;

