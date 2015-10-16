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


function Data_Analysis_List2_9_5(nFile)

    if nargin<1
        nFile = 1;
    end
        
    addpath('../Func');
    setDir;    
    fileDirName         = fileDirNames{nFile}; %#ok<NASGU,USENS>
    fileName            = fileNames{nFile}; %#ok<USENS>
    load([tempDatDir, fileName, '.mat']);
    load([tempDatDir, fileName, '_HalfTimeBilateralIndex.mat']);
    
    RSquareThres        = 0.5;
    nTime               = length(timePoints)-1;
    xTrack              = squeeze(mean(tracks(:, timePoints(nTime)+1:timePoints(nTime+1), 1), 2)); %#ok<NODEF>
    yTrack              = squeeze(mean(tracks(:, timePoints(nTime)+1:timePoints(nTime+1), 2), 2));
    
    figure;
    hold on;
    plot(xTrack(RSquare<RSquareThres), yTrack(RSquare<RSquareThres), 'o', 'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5])
    scatter(xTrack(RSquare>=RSquareThres), yTrack(RSquare>=RSquareThres), [], halfTime(RSquare>=RSquareThres), 'fill')
    hold off;
    ylim([000 400])
    xlim([0 1600])
    box off
    xlabel('Last time R-C location (um)')
    ylabel('Last time M-L location (um)')
    colorbar
    setPrint(8, 6, [plotDir, 'BilateralIndexHalfTimeThre05_', fileName], 'pdf')    
end