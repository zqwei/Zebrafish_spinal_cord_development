%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.  Covariance analysis: clustering analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

function Data_Analysis_List1_1n(nFile, nTimes)
    % load data
    addpath('../Func');
    setDir;
    
    fileDirName       = fileDirNames{nFile};
    fileName          = fileNames{nFile};
    load([tempDatDir, fileName, '.mat']);
    %
    %

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1.1  Covariance analysis: plotting the covariance matrix every T
    % timepoints by removing the points with weak linear correlations with the
    % other points in the neuronal pool
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %
    dff           = dffWithNoDuplicates;
%     ind           = indWithNoDuplicates;
    tracks        = tracksWithNoDuplicates;
    idx           = idxWithNoDuplicates;
    clear dffWithNoDuplicates indWithNoDuplicates;
    timeStep      = 1200;
    numT          = size(dff, 2);
    numPlot       = ceil(numT/timeStep);
    timePoints    = [0 (1:(numPlot-1))*timeStep numT];
    m             = length(nTimes);
    corrThres     = 0.4; 



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1.1.2  Covariance analysis: plotting the covariance matrix every T
    % timepoints by removing the points with weak linear correlations with the
    % other points in the neuronal pool -- using local leaf order at last time
    % point
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure;    
    nPlot           = numPlot;
    slicedDFF       = dff(:,timePoints(nPlot-9)+1:timePoints(nPlot+1));  
    % The length of average also needs to be check      
    slicedCorr      = corr(slicedDFF');
    tSlicedCorr     = slicedCorr - eye(size(slicedCorr));
    slicedIndex     = max(tSlicedCorr,[],2)>corrThres;% max(abs(tSlicedCorr),[],2)>corrThres;
    slicedDFF       = slicedDFF(slicedIndex, :);
    distNeurons     = pdist(slicedDFF, 'correlation');
    linkNeurons     = linkage(slicedDFF,'single','correlation');
    leafOrder       = optimalleaforder(linkNeurons, distNeurons);       

    for nTime           = 1:m
        subplot(1, m, nTime)
        nPlot           = nTimes(nTime);
        slicedDFF       = dff(slicedIndex,timePoints(nPlot)+1:timePoints(nPlot+1));        
        imagesc(corr(slicedDFF(leafOrder,:)'),[-1 1]);
        xlabel('Neuronal index')
        ylabel('Neuronal index')
        title(['Time from ' num2str((nPlot-1)*timeStep/1000) 'k to ' num2str(nPlot*timeStep/1000) 'k']);
        box off;
    end
    setFigureSize(m*8, 6);
    
    
    
    
    
    figure;    
    nPlot           = numPlot;
    slicedDFF       = dff(:,timePoints(nPlot-9)+1:timePoints(nPlot+1));  
    % The length of average also needs to be check      
    slicedCorr      = corr(slicedDFF');
    tSlicedCorr     = slicedCorr - eye(size(slicedCorr));
    slicedIndex     = max(tSlicedCorr,[],2)>corrThres;% max(abs(tSlicedCorr),[],2)>corrThres;
    for nTime           = 1:m
        subplot(1, m, nTime)
        nPlot           = nTimes(nTime);
        slicedDFF       = dff(slicedIndex,timePoints(nPlot)+1:timePoints(nPlot+1));  
        distNeurons     = pdist(slicedDFF, 'correlation');
        linkNeurons     = linkage(slicedDFF,'single','correlation');
        leafOrder       = optimalleaforder(linkNeurons, distNeurons);       
        imagesc(corr(slicedDFF(leafOrder,:)'),[-1 1]);
        xlabel('Neuronal index')
        ylabel('Neuronal index')
        title(['Time from ' num2str((nPlot-1)*timeStep/1000) 'k to ' num2str(nPlot*timeStep/1000) 'k']);
        box off;
    end
    setFigureSize(m*8, 6);


end