%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.  Explained variance: Factor analysis -- corrected loading matrix --
% sigmoid fit
% 
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function FAEV_v0_1(nFile)            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>    
    mnx               = [];
    load([tempDatDir, fileName, '.mat'], 'mnx', 'tracks', 'side', 'leafOrder', 'slicedIndex'); 
    load([tempDatDir, 'EV_' fileName, '.mat'], 'EVMat')
    
    numTime           = size(EVMat, 2); %#ok<USENS>
    numNeuron         = size(EVMat, 1);
    meantracks        = squeeze(mean(tracks, 2)); %#ok<NODEF>
    halfEVTime        = nan(numNeuron, 1);
    RSquare           = zeros(numNeuron, 1);
    
    m                 = ceil(sqrt(numNeuron));
    
    neuronName        = find(slicedIndex);
    neuronName        = neuronName(leafOrder);
    
    figure;
    
    for nNeuron       = 1:numNeuron
        subplot(m, m, nNeuron)
        
        if side(nNeuron) == 1
            mColor  = [0    0.4470    0.7410];
        else
            mColor  = [0.6350    0.0780    0.1840];
        end
        
        timeIndex = ~isnan(EVMat(nNeuron, :));
        if sum(timeIndex) > 10
            [fitParams, fitResult] = sigm_fit(find(timeIndex)/60, EVMat(nNeuron, timeIndex)');
            RSquare(nNeuron)         = 1 - mean((EVMat(nNeuron, timeIndex)' - fitResult.ypred).^2)./var(EVMat(nNeuron, timeIndex)');
            halfEVTime(nNeuron)      = fitParams(3);
        end
        
        hold on;
        plot(find(timeIndex)/60, EVMat(nNeuron, timeIndex), 'o', 'color', mColor);
        if ~isnan(halfEVTime(nNeuron))
            plot(find(timeIndex)/60, fitResult.ypred, '-', 'linewid', 2.0, 'Color', mColor);
            plot(find(timeIndex)/60, fitResult.ypredlowerCI, '-', 'linewid', 0.5, 'Color', mColor); 
            plot(find(timeIndex)/60, fitResult.ypredupperCI, '-', 'linewid', 0.5, 'Color', mColor); 
        end
        hold off;
        
        title(['#' num2str(neuronName(nNeuron)) ' R^2=' num2str(RSquare(nNeuron))])
        
        if ~isempty(mnx)
            if mnx(nNeuron)
                title(['#' num2str(neuronName(nNeuron)) ' MNX+; R^2=' num2str(RSquare(nNeuron))])
            else
                title(['#' num2str(neuronName(nNeuron)) ' MNX-; R^2=' num2str(RSquare(nNeuron))])
            end
        end
        
        ylim([0 1])
        xlim([0 numTime/60])    
        xlabel('Time (hour)')
        ylabel('EV')
        box off
    end
    
    setPrint(8*m, 6*m, [plotDir, 'SigFitEV_', fileName], 'pdf');
    
    RSquareThres = 0.9;
    
    figure;
    xTrack       = meantracks(:, 1);
    yTrack       = meantracks(:, 2);
    zTrack       = meantracks(:, 3);
    
    halfEVTime(halfEVTime<0) = 0;
    figure;
    hold on
    plot(xTrack(RSquare<RSquareThres), yTrack(RSquare<RSquareThres), 'o', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', [0.5 0.5 0.5] );
%     scatter(xTrack(RSquare>=RSquareThres), yTrack(RSquare>=RSquareThres), [], halfEVTime(RSquare>=RSquareThres), 'fill')
    
    if isempty(mnx)
        scatter(xTrack(RSquare>=RSquareThres), yTrack(RSquare>=RSquareThres), [], halfEVTime(RSquare>=RSquareThres), 'fill')
    else
        if sum(RSquare>=RSquareThres & mnx)>0
            scatter(xTrack(RSquare>=RSquareThres & mnx), yTrack(RSquare>=RSquareThres & mnx), [], halfEVTime(RSquare>=RSquareThres & mnx), 'o', 'fill')
        end
        hold on
        if sum(RSquare>=RSquareThres & ~mnx)>0
            scatter(xTrack(RSquare>=RSquareThres & ~mnx), yTrack(RSquare>=RSquareThres & ~mnx), [], halfEVTime(RSquare>=RSquareThres & ~mnx), 's', 'fill')
        end
        hold off
    end
    hold off
    
    xlabel('x-axis')
    ylabel('y-axis')
    box off
    colorbar
        
    setPrint(8, 6, [plotDir, 'HalfEVTimeXY_', fileName], 'pdf');
    
    figure;
    hold on
    plot(xTrack(RSquare<RSquareThres), zTrack(RSquare<RSquareThres), 'o', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', [0.5 0.5 0.5] );
%     scatter(xTrack(RSquare>=RSquareThres), zTrack(RSquare>=RSquareThres), [], halfEVTime(RSquare>=RSquareThres), 'fill')
    
    if isempty(mnx)
        scatter(xTrack(RSquare>=RSquareThres), zTrack(RSquare>=RSquareThres), [], halfEVTime(RSquare>=RSquareThres), 'fill')
    else
        if sum(RSquare>=RSquareThres & mnx)>0
            scatter(xTrack(RSquare>=RSquareThres & mnx), zTrack(RSquare>=RSquareThres & mnx), [], halfEVTime(RSquare>=RSquareThres & mnx), 'o', 'fill')
        end
        hold on
        if sum(RSquare>=RSquareThres & ~mnx)>0
            scatter(xTrack(RSquare>=RSquareThres & ~mnx), zTrack(RSquare>=RSquareThres & ~mnx), [], halfEVTime(RSquare>=RSquareThres & ~mnx), 's', 'fill')
        end
        hold off
    end
    hold off
    
    xlabel('x-axis')
    ylabel('z-axis')
    box off
    colorbar
        
    setPrint(8, 6, [plotDir, 'HalfEVTimeXZ_', fileName], 'pdf');
    
    close all
    save([tempDatDir, 'EV_' fileName, '.mat'], 'halfEVTime', 'RSquare', '-append');
end