%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.  Covariance analysis: clustering analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

function Data_Analysis_List1_1(nFile)
    % load data
    addpath('../Func');
    setDir;
    
    fileDirName       = fileDirNames{nFile};
    fileName          = fileNames{nFile};
    % Data names:
    % dffWithNoDuplicates
    % indWithNoDuplicates
    % tracksWithNoDuplicates
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
    m             = ceil(sqrt(numPlot + 1));
    corrThres     = 0.4; 
    % Here corrThres is chosen randomly, one should test if all neurons below 
    % this corrThres is mostly the non-selective neurons as those predefined

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1.1.1  Covariance analysis: plotting the covariance matrix every T
    % timepoints by removing the points with weak linear correlations with the
    % other points in the neuronal pool -- using local leaf order
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % nskip             = true;
    % if ~nskip
    %     nGroup        = 2;
    %     minSizeGroup  = 3;    
    %     figure;
    %     for nPlot           = 1:numPlot
    %         subplot(m, m, nPlot)
    %         slicedDFF       = dff(:,timePoints(nPlot)+1:timePoints(nPlot+1));        
    %         slicedCorr      = corr(slicedDFF');
    %         tSlicedCorr     = slicedCorr - eye(size(slicedCorr));
    %         slicedIndex     = max(tSlicedCorr,[],2)>corrThres;
    %         slicedDFF       = slicedDFF(slicedIndex, :);
    %         distNeurons     = pdist(slicedDFF, 'correlation');
    %         linkNeurons     = linkage(slicedDFF,'single','correlation');
    %         leafOrder       = optimalleaforder(linkNeurons, distNeurons);
    %         imagesc(corr(slicedDFF(leafOrder,:)'),[-1 1]);
    %         xlabel('Neuronal index')
    %         ylabel('Neuronal index')
    %         title(['Time from ' num2str((nPlot-1)*timeStep/1000) 'k to ' num2str(nPlot*timeStep/1000) 'k']);
    %         box off;
    %     end
    %     setPrint(m*8, m*6, [plotDir, 'KnockingOut_UngroupedNeurons_', fileName], 'pdf');
    % end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1.1.2  Covariance analysis: plotting the covariance matrix every T
    % timepoints by removing the points with weak linear correlations with the
    % other points in the neuronal pool -- using local leaf order at last time
    % point
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nskip             = false;
    if ~nskip
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%%% 3d information of weak corrlated cells
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     track_x         = squeeze(mean(tracks(:, timePoints(nPlot)+1:timePoints(nPlot+1),1), 2));
    %     track_y         = squeeze(mean(tracks(:, timePoints(nPlot)+1:timePoints(nPlot+1),2), 2));
    %     track_z         = squeeze(mean(tracks(:, timePoints(nPlot)+1:timePoints(nPlot+1),3), 2));
    %     var_x           = squeeze(var(tracks(:, timePoints(nPlot)+1:timePoints(nPlot+1),1), [], 2));
    %     var_y           = squeeze(var(tracks(:, timePoints(nPlot)+1:timePoints(nPlot+1),2), [], 2));
    %     var_z           = squeeze(var(tracks(:, timePoints(nPlot)+1:timePoints(nPlot+1),3), [], 2));
    %     radius          = sqrt(var_x + var_y + var_z);
    %     coordinates                = squeeze(mean(tracks(:, timePoints(nPlot)+1:timePoints(nPlot+1),:), 2));
    %     colorTable                 = zeros(size(coordinates));
    %     colorTable(slicedIndex,1)  = 1;
    %     colorTable(~slicedIndex,3) = 1;
    %     
    %     flist                      = dir([fileDirName '/image_data/Period_TM*.tif']);
    % 
    %     exportMask([fileDirName '/image_data/' flist(end).name], [plotDir 'Overlay_' fileName], coordinates, colorTable)
        for nPlot           = 1:numPlot
            subplot(m, m, nPlot)
            slicedDFF       = dff(slicedIndex,timePoints(nPlot)+1:timePoints(nPlot+1));        
            imagesc(corr(slicedDFF(leafOrder,:)'),[-1 1]);
            xlabel('Neuronal index')
            ylabel('Neuronal index')
            title(['Time from ' num2str((nPlot-1)*timeStep/1000) 'k to ' num2str(nPlot*timeStep/1000) 'k']);
            box off;
        end
        setPrint(m*8, m*6, [plotDir, 'CorrMat_', fileName], 'pdf');
    end
    
    
    
    
    
    nskip             = false;
    if ~nskip
        figure;    
        nPlot           = numPlot;
        slicedDFF       = dff(:,timePoints(nPlot-9)+1:timePoints(nPlot+1));  
        % The length of average also needs to be check      
        slicedCorr      = corr(slicedDFF');
        tSlicedCorr     = slicedCorr - eye(size(slicedCorr));
        slicedIndex     = max(tSlicedCorr,[],2)>corrThres;% max(abs(tSlicedCorr),[],2)>corrThres;
        for nPlot           = 1:numPlot
            subplot(m, m, nPlot)
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
        setPrint(m*8, m*6, [plotDir, 'CorrMatLocal_', fileName], 'pdf');
    end

%     nPlot           = numPlot;
%     slicedDFF       = dff(:,timePoints(nPlot-9)+1:timePoints(nPlot+1));  
%     % The length of average also needs to be check      
%     slicedCorr      = corr(slicedDFF');
%     tSlicedCorr     = slicedCorr - eye(size(slicedCorr));
%     slicedIndex     = max(tSlicedCorr,[],2)>corrThres;% max(abs(tSlicedCorr),[],2)>corrThres;
%     slicedDFF       = slicedDFF(slicedIndex, :);
%     distNeurons     = pdist(slicedDFF, 'correlation');
%     linkNeurons     = linkage(slicedDFF,'single','correlation');
%     leafOrder       = optimalleaforder(linkNeurons, distNeurons);
% 
%     dff           = dff(slicedIndex,:);
%     tracks        = tracks(slicedIndex,:,:);
%     idx           = idx(slicedIndex,:);
% 
%     dff           = dff(leafOrder,:);
%     tracks        = tracks(leafOrder,:,:);
%     idx           = idx(leafOrder,:);
% 
    save([tempDatDir, fileName, '.mat'], 'dff', 'tracks','idx', 'leafOrder', 'slicedIndex','-append');
end