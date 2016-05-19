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


function FAEV_v0_3(nFile)            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>    
    mnx               = [];
    load([tempDatDir, fileName, '.mat'], 'mnx', 'tracks', 'side', 'leafOrder', 'slicedIndex'); 
    load([tempDatDir, 'EV_' fileName, '.mat'], 'EVMat', 'halfEVTime', 'RSquare', 'halfActTime', 'validFitIndex')
    
    meantracks        = squeeze(mean(tracks, 2));         
    RSquareThres      = 0.6;
    
    figure;
    xTrack       = meantracks(:, 1);
    yTrack       = meantracks(:, 2);
    zTrack       = meantracks(:, 3);
    
    figure;
    hold on
    plot(xTrack(RSquare<RSquareThres | ~validFitIndex), yTrack(RSquare<RSquareThres | ~validFitIndex), 'o', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', [0.5 0.5 0.5] );
%     scatter(xTrack(RSquare>=RSquareThres), yTrack(RSquare>=RSquareThres), [], halfEVTime(RSquare>=RSquareThres), 'fill')
    
    if isempty(mnx)
        scatter(xTrack(RSquare>=RSquareThres & validFitIndex), yTrack(RSquare>=RSquareThres & validFitIndex), [], halfEVTime(RSquare>=RSquareThres & validFitIndex), 'fill')
    else
        if sum(RSquare>=RSquareThres & mnx & validFitIndex)>0
            scatter(xTrack(RSquare>=RSquareThres & mnx & validFitIndex), yTrack(RSquare>=RSquareThres & mnx & validFitIndex), [], halfEVTime(RSquare>=RSquareThres & mnx & validFitIndex), 'o', 'fill')
        end
        hold on
        if sum(RSquare>=RSquareThres & ~mnx & validFitIndex)>0
            scatter(xTrack(RSquare>=RSquareThres & ~mnx & validFitIndex), yTrack(RSquare>=RSquareThres & ~mnx & validFitIndex), [], halfEVTime(RSquare>=RSquareThres & ~mnx & validFitIndex), 's', 'fill')
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
    plot(xTrack(RSquare<RSquareThres | ~validFitIndex), zTrack(RSquare<RSquareThres | ~validFitIndex), 'o', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', [0.5 0.5 0.5] );    
    if isempty(mnx)
        scatter(xTrack(RSquare>=RSquareThres & validFitIndex), zTrack(RSquare>=RSquareThres & validFitIndex), [], halfEVTime(RSquare>=RSquareThres & validFitIndex), 'fill')
    else
        if sum(RSquare>=RSquareThres & mnx & validFitIndex)>0
            scatter(xTrack(RSquare>=RSquareThres & mnx & validFitIndex), zTrack(RSquare>=RSquareThres & mnx & validFitIndex), [], halfEVTime(RSquare>=RSquareThres & mnx & validFitIndex), 'o', 'fill')
        end
        hold on
        if sum(RSquare>=RSquareThres & ~mnx & validFitIndex)>0
            scatter(xTrack(RSquare>=RSquareThres & ~mnx & validFitIndex), zTrack(RSquare>=RSquareThres & ~mnx & validFitIndex), [], halfEVTime(RSquare>=RSquareThres & ~mnx & validFitIndex), 's', 'fill')
        end
        hold off
    end
    hold off
    
    xlabel('x-axis')
    ylabel('z-axis')
    box off
    colorbar
        
    setPrint(8, 6, [plotDir, 'HalfEVTimeXZ_', fileName], 'pdf');
    
    
    
    figure;
    hold on
    plot(xTrack(RSquare<RSquareThres | ~validFitIndex), yTrack(RSquare<RSquareThres | ~validFitIndex), 'o', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', [0.5 0.5 0.5] );
%     scatter(xTrack(RSquare>=RSquareThres), yTrack(RSquare>=RSquareThres), [], halfEVTime(RSquare>=RSquareThres), 'fill')
    
    if isempty(mnx)
        scatter(xTrack(RSquare>=RSquareThres & validFitIndex), yTrack(RSquare>=RSquareThres & validFitIndex), [], halfActTime(RSquare>=RSquareThres & validFitIndex), 'fill')
    else
        if sum(RSquare>=RSquareThres & mnx & validFitIndex)>0
            scatter(xTrack(RSquare>=RSquareThres & mnx & validFitIndex), yTrack(RSquare>=RSquareThres & mnx & validFitIndex), [], halfActTime(RSquare>=RSquareThres & mnx & validFitIndex), 'o', 'fill')
        end
        hold on
        if sum(RSquare>=RSquareThres & ~mnx & validFitIndex)>0
            scatter(xTrack(RSquare>=RSquareThres & ~mnx & validFitIndex), yTrack(RSquare>=RSquareThres & ~mnx & validFitIndex), [], halfActTime(RSquare>=RSquareThres & ~mnx & validFitIndex), 's', 'fill')
        end
        hold off
    end
    hold off
    
    xlabel('x-axis')
    ylabel('y-axis')
    box off
    colorbar
        
    setPrint(8, 6, [plotDir, 'HalfActTimeXY_', fileName], 'pdf');
    
    figure;
    hold on
    plot(xTrack(RSquare<RSquareThres | ~validFitIndex), zTrack(RSquare<RSquareThres | ~validFitIndex), 'o', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', [0.5 0.5 0.5] );    
    if isempty(mnx)
        scatter(xTrack(RSquare>=RSquareThres & validFitIndex), zTrack(RSquare>=RSquareThres & validFitIndex), [], halfActTime(RSquare>=RSquareThres & validFitIndex), 'fill')
    else
        if sum(RSquare>=RSquareThres & mnx & validFitIndex)>0
            scatter(xTrack(RSquare>=RSquareThres & mnx & validFitIndex), zTrack(RSquare>=RSquareThres & mnx & validFitIndex), [], halfActTime(RSquare>=RSquareThres & mnx & validFitIndex), 'o', 'fill')
        end
        hold on
        if sum(RSquare>=RSquareThres & ~mnx & validFitIndex)>0
            scatter(xTrack(RSquare>=RSquareThres & ~mnx & validFitIndex), zTrack(RSquare>=RSquareThres & ~mnx & validFitIndex), [], halfActTime(RSquare>=RSquareThres & ~mnx & validFitIndex), 's', 'fill')
        end
        hold off
    end
    hold off
    
    xlabel('x-axis')
    ylabel('z-axis')
    box off
    colorbar
        
    setPrint(8, 6, [plotDir, 'HalfActTimeXZ_', fileName], 'pdf');

end