%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  Covariance analysis: factor analysis
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.4 FA vs islands -- clustering
%     following 2.7
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function Data_Analysis_List2_8_4(nFile)

    if nargin<1
        nFile = 1;
    end
    
    pixelToImageRatio = 0.406;
    scaleBarOffSet    = 5;
    scaleBarLength    = 50;
    
    addpath('../Func');
    setDir;    
    fileDirName         = fileDirNames{nFile}; %#ok<NASGU,USENS>
    fileName            = fileNames{nFile}; %#ok<USENS>
    load([tempDatDir, fileName, '.mat']);
    load([tempDatDir, fileName, '_HalfTime.mat']);
    
    RSquareThres        = 0.6;
    % EVThres             = 0.5; % maybe also consider EV above 0.5
    nTime               = length(timePoints)-1;
    xTrack              = squeeze(mean(tracks(:, timePoints(nTime)+1:timePoints(nTime+1), 1), 2)); %#ok<NODEF>
    yTrack              = squeeze(mean(tracks(:, timePoints(nTime)+1:timePoints(nTime+1), 2), 2));
    zTrack              = squeeze(mean(tracks(:, timePoints(nTime)+1:timePoints(nTime+1), 3), 2)); %#ok<NASGU>
    
    figure;
    hold on;
    plot(xTrack(RSquare<RSquareThres | idx == 3), yTrack(RSquare<RSquareThres | idx == 3), 'o', 'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5])
    scatter(xTrack(RSquare>=RSquareThres & idx ~= 3), yTrack(RSquare>=RSquareThres  & idx ~= 3), [], halfTime(RSquare>=RSquareThres & idx ~= 3), 'fill')
    
    plot([scaleBarOffSet scaleBarOffSet+scaleBarLength/pixelToImageRatio], [scaleBarOffSet scaleBarOffSet], '-k', 'linewid', 2.0)
    plot([scaleBarOffSet scaleBarOffSet], [scaleBarOffSet+scaleBarLength/pixelToImageRatio scaleBarOffSet], '-k', 'linewid', 2.0)
    
    hold off;
    ylim([000 400])
    xlim([0 1600])
%     caxis([0 3.5])
    axis off
    box off
    xlabel('Last time R-C location (um)')
    ylabel('Last time M-L location (um)')
    colorbar
    setPrint(8, 6, [plotDir, 'HalfTimeThre05_', fileName], 'pdf')  
    
    
%     figure;
%     hold on;
% %     plot3(xTrack(RSquare<RSquareThres | idx == 3),...
% %           yTrack(RSquare<RSquareThres | idx == 3),...
% %           zTrack(RSquare<RSquareThres | idx == 3)*6, ...
% %           'o', 'MarkerEdgeColor', [0.5 0.5 0.5],...
% %           'MarkerFaceColor', [0.5 0.5 0.5])
% %     scatter3(xTrack(RSquare>=RSquareThres & idx ~= 3), yTrack(RSquare>=RSquareThres  & idx ~= 3), zTrack(RSquare>=RSquareThres  & idx ~= 3)*6, [], halfTime(RSquare>=RSquareThres), 'fill')
%     scatter3(xTrack(RSquare>=RSquareThres & idx ~= 3 & mnx(orgIdx) == 0), ...
%           yTrack(RSquare>=RSquareThres  & idx ~= 3 & mnx(orgIdx) == 0),...
%           zTrack(RSquare>=RSquareThres  & idx ~= 3 & mnx(orgIdx) == 0)*6,...
%           [], halfTime(RSquare>=RSquareThres  & idx ~= 3 & mnx(orgIdx) == 0))
%       
%     scatter3(xTrack(RSquare>=RSquareThres & idx ~= 3 & mnx(orgIdx) == 1), ...
%           yTrack(RSquare>=RSquareThres  & idx ~= 3 & mnx(orgIdx) == 1),...
%           zTrack(RSquare>=RSquareThres  & idx ~= 3 & mnx(orgIdx) == 1)*6,...
%           [], halfTime(RSquare>=RSquareThres  & idx ~= 3 & mnx(orgIdx) == 1), 'fill')  
%     hold off;
%     ylim([000 400])
%     xlim([000 1600])
%     box off
%     xlabel('Last time R-C location (um)')
%     ylabel('Last time M-L location (um)')
%     colorbar
%     setPrint(8, 6, [plotDir, 'HalfTimeThre05_', fileName], 'pdf')
    
    
    
%     figure;
%     nTimeStep            = length(timePoints)-1;
%     nPlot                = 0;
%     for nTime            = nTimeStep-14:7:nTimeStep
%         nPlot            = nPlot + 1;
%         subplot(1, 3, nPlot)
%         plot((timePoints(nTime)+1:timePoints(nTime+1))/4/3600, dff(44,timePoints(nTime)+1:timePoints(nTime+1)));
%         xlim([timePoints(nTime)+1 timePoints(nTime+1)]/4/3600)
%         xlabel('Time (hour)')
%         ylabel('df/f')
%         box off
%     end
%     setPrint(8*3, 6, [plotDir, 'ExampleNeuronIndex44'], 'pdf')
    
    
    
    
end