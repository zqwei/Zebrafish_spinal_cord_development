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


function Data_Analysis_List2_9_6(nFile)

    if nargin<1
        nFile = 1;
    end
        
    addpath('../Func');
    setDir;    
    fileDirName         = fileDirNames{nFile}; %#ok<NASGU,USENS>
    fileName            = fileNames{nFile}; %#ok<USENS>
    load([tempDatDir, fileName, '.mat']);
    load([tempDatDir, fileName, '_HalfTimeBilateralIndex.mat']);
    load([tempDatDir, fileName, '_FAEV.mat'], 'opt2Dim');

    RSquareThres        = 0.6;        
    BIHalfTime          = abs(halfTime);
    BIplotNeurons       = RSquare>=RSquareThres & BIHalfTime>timePoints(1)/4/3600 & BIHalfTime<timePoints(end-1)/4/3600; %#ok<COLND>
        
    load([tempDatDir, fileName, '_HalfTime.mat']);
    EVHalfTime          = halfTime;
    EVplotNeurons       = RSquare>=RSquareThres & EVHalfTime>timePoints(1)/4/3600 & EVHalfTime<timePoints(end-1)/4/3600; %#ok<COLND>    
    neuronToPlot        = find(BIplotNeurons & EVplotNeurons);
    
    
    
    for nPlot           = 1:length(neuronToPlot);
        figure;
        nNeuron         = neuronToPlot(nPlot);        
        subplot(2, 2, 1 );
        nTime           = sum(timePoints/4/3600<EVHalfTime(nNeuron));
        refDFF          = dff(:,timePoints(nTime)+1:timePoints(nTime+1)); %#ok<NODEF>
        refDFF          = bsxfun(@minus, refDFF, mean(refDFF,2));
        refDFF          = bsxfun(@rdivide, refDFF, std(refDFF,[],2))';
        [lambda,psi]    = factoran(refDFF, opt2Dim(nTime), 'scores','regression', 'rotate', 'none');
        DFFEst          = LONOFASingleUnit (refDFF, lambda, psi, nNeuron);
        hold on;
        plot((timePoints(nTime)+1:timePoints(nTime+1))/4/3600,refDFF(:,nNeuron),'-b','linewid',0.5);
        plot((timePoints(nTime)+1:timePoints(nTime+1))/4/3600,DFFEst,'-r','linewid',0.5);
        hold off
        xlim([(timePoints(nTime)+1)/4/3600 timePoints(nTime+1)/4/3600])
        xlabel('Time (hour)')
        ylabel('Normalized DF/F')
        title('EV half time')
        box off

        subplot(2, 2, 2);
        nTime           = sum(timePoints/4/3600<BIHalfTime(nNeuron));
        refDFF          = dff(:,timePoints(nTime)+1:timePoints(nTime+1));
        refDFF          = bsxfun(@minus, refDFF, mean(refDFF,2));
        refDFF          = bsxfun(@rdivide, refDFF, std(refDFF,[],2))';
        [lambda,psi]    = factoran(refDFF, opt2Dim(nTime), 'scores','regression', 'rotate', 'none');
        DFFEst          = LONOFASingleUnit (refDFF, lambda, psi, nNeuron);
        hold on;
        plot((timePoints(nTime)+1:timePoints(nTime+1))/4/3600,refDFF(:,nNeuron),'-b','linewid',0.5);
        plot((timePoints(nTime)+1:timePoints(nTime+1))/4/3600,DFFEst,'-r','linewid',0.5);
        hold off
        xlim([(timePoints(nTime)+1)/4/3600 timePoints(nTime+1)/4/3600])
        xlabel('Time (hour)')
        ylabel('Normalized DF/F')
        title('BI saturation time')
        box off
        
        
        subplot(2, 2, 3 );
        nTime           = sum(timePoints/4/3600<EVHalfTime(nNeuron));
        SlicedDFF       = dff(:,timePoints(nTime)+1:timePoints(nTime+1));
        corrDFF         = corr(SlicedDFF');
        corrDFFnNeuron  = corrDFF(nNeuron,:);
        xTrack          = squeeze(mean(tracks(:, timePoints(nTime)+1:timePoints(nTime+1), 1), 2)); %#ok<NODEF>
        yTrack          = squeeze(mean(tracks(:, timePoints(nTime)+1:timePoints(nTime+1), 2), 2));
        scatter(xTrack, yTrack, [], corrDFFnNeuron, '.')
        colorbar
        caxis([-1 1])
        xlim([(timePoints(nTime)+1)/4/3600 timePoints(nTime+1)/4/3600])
        ylim([000 400])
        xlim([0 1600])
        box off
        xlabel('Last time R-C location (um)')
        ylabel('Last time M-L location (um)')
        title('EV half time')
        
        subplot(2, 2, 4 );
        nTime           = sum(timePoints/4/3600<BIHalfTime(nNeuron));
        SlicedDFF       = dff(:,timePoints(nTime)+1:timePoints(nTime+1));
        corrDFF         = corr(SlicedDFF');
        corrDFFnNeuron  = corrDFF(nNeuron,:);
        xTrack          = squeeze(mean(tracks(:, timePoints(nTime)+1:timePoints(nTime+1), 1), 2));
        yTrack          = squeeze(mean(tracks(:, timePoints(nTime)+1:timePoints(nTime+1), 2), 2));
        scatter(xTrack, yTrack, [], corrDFFnNeuron, '.')
        colorbar
        caxis([-1 1])
        xlim([(timePoints(nTime)+1)/4/3600 timePoints(nTime+1)/4/3600])
        ylim([000 400])
        xlim([0 1600])
        box off
        xlabel('Last time R-C location (um)')
        ylabel('Last time M-L location (um)')
        title('BI saturation time')        
        
        setPrint(8*2, 6*2, [plotDir, 'ExampleNeuron', num2str(nNeuron,'%02d'),'_', fileName], 'pdf')            
        close all
    end
end