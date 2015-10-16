%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.  Covariance analysis: cross-middle line activity neurons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 3.1 Analysis the correlation matrix based on the first few time point
% order
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
%

function Data_Analysis_List3_1(nFile)
    % load data
    addpath('../Func');
    setDir;
    
    fileDirName       = fileDirNames{nFile};
    fileName          = fileNames{nFile};

    load([tempDatDir, fileName, '.mat']);

    dff           = dffWithNoDuplicates;
    clear dffWithNoDuplicates;
    timeStep      = 1200;
    numT          = size(dff, 2);
    numPlot       = ceil(numT/timeStep);
    timePoints    = [0 (1:(numPlot-1))*timeStep numT];
    m             = ceil(sqrt(numPlot + 1));
    corrThres     = 0.4; 

    figure;    
    nPlot           = numPlot;
    slicedDFF       = dff(:,timePoints(nPlot-9)+1:timePoints(nPlot+1));  
    % The length of average also needs to be check      
    slicedCorr      = corr(slicedDFF');
    tSlicedCorr     = slicedCorr - eye(size(slicedCorr));
    slicedIndex     = max(tSlicedCorr,[],2)>corrThres;
    
    slicedDFF       = dff(:,timePoints(2)+1:timePoints(6));
    slicedDFF       = slicedDFF(slicedIndex, :);
    distNeurons     = pdist(slicedDFF, 'correlation');
    linkNeurons     = linkage(slicedDFF,'single','correlation');
    leafOrder       = optimalleaforder(linkNeurons, distNeurons);
    
    for nPlot           = 1:numPlot
        subplot(m, m, nPlot)
        slicedDFF       = dff(slicedIndex,timePoints(nPlot)+1:timePoints(nPlot+1));  
        imagesc(corr(slicedDFF(leafOrder,:)'),[-1 1]);
        xlabel('Neuronal index')
        ylabel('Neuronal index')
        title(['Time from ' num2str((nPlot-1)*timeStep/1000) 'k to ' num2str(nPlot*timeStep/1000) 'k']);
        box off;
    end
    setPrint(m*8, m*6, [plotDir, 'CorrMatFirst_', fileName], 'pdf');

end