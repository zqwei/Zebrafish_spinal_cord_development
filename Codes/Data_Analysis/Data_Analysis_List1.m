%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.  Covariance analysis: clustering analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

function Data_Analysis_List1(nFile)
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
    % timepoints
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %
    dff           = dffWithNoDuplicates;
%     ind           = indWithNoDuplicates;
    tracks        = tracksWithNoDuplicates;
    clear dffWithNoDuplicates indWithNoDuplicates;

    timeStep      = 1200;
    numT          = size(dff, 2);
    numPlot       = ceil(numT/timeStep);
    timePoints    = [0 (1:(numPlot-1))*timeStep numT];

    m             = ceil(sqrt(numPlot + 1));


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1.1.1 plot with orginial clustering
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nskip             = true;
    if ~nskip
        figure;
        for nPlot     = 1:numPlot
            subplot(m, m, nPlot)
            slicedDFF = dff(:,timePoints(nPlot)+1:timePoints(nPlot+1));
            imagesc(corr(slicedDFF'),[-1 1]);
            xlabel('Neuronal index')
            ylabel('Neuronal index')
            title(['Time from ' num2str((nPlot-1)*timeStep/1000) 'k to ' num2str(nPlot*timeStep/1000) 'k']);
            box off;
        end
        % across time points
        subplot(m, m, numPlot + 1)
        slicedDFF = dff;
        imagesc(corr(slicedDFF'),[-1 1]);
        xlabel('Neuronal index')
        ylabel('Neuronal index')
        title(['Time from ' num2str(0*timeStep/1000) 'k to ' num2str(nPlot*timeStep/1000) 'k']);
        box off;
        setPrint(m*8, m*6, [plotDir, 'Correlation_OrginialClustering_', fileName], 'pdf');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1.1.2 plot with reindexing from clustering
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nskip             = true;
    if ~nskip
        figure;
        slicedDFF           = dff(:,timePoints(numPlot)+1:timePoints(numPlot+1));
        distNeurons         = pdist(slicedDFF, 'correlation');
        linkNeurons         = linkage(slicedDFF,'single','correlation');
        leafOrder           = optimalleaforder(linkNeurons, distNeurons);
        for nPlot           = 1:numPlot
            subplot(m, m, nPlot)
            slicedDFF       = dff(:,timePoints(nPlot)+1:timePoints(nPlot+1));
    %         distNeurons     = pdist(slicedDFF, 'correlation');
    %         linkNeurons     = linkage(slicedDFF,'single','correlation');
    %         leafOrder       = optimalleaforder(linkNeurons, distNeurons);
            imagesc(corr(slicedDFF(leafOrder,:)'),[-1 1]);
            xlabel('Neuronal index')
            ylabel('Neuronal index')
            title(['Time from ' num2str((nPlot-1)*timeStep/1000) 'k to ' num2str(nPlot*timeStep/1000) 'k']);
        end
        % across time points
        subplot(m, m, numPlot + 1)
        slicedDFF = dff;
    %     distNeurons     = pdist(slicedDFF, 'correlation');
    %     linkNeurons     = linkage(slicedDFF,'single','correlation');
    %     leafOrder       = optimalleaforder(linkNeurons, distNeurons);
        imagesc(corr(slicedDFF(leafOrder,:)'),[-1 1]);
        xlabel('Neuronal index')
        ylabel('Neuronal index')
        title(['Time from ' num2str(0*timeStep/1000) 'k to ' num2str(nPlot*timeStep/1000) 'k']);
        box off;
        setPrint(m*8, m*6, [plotDir, 'Correlation_ModifiedClustering_', fileName], 'pdf');
    end



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1.2  Covariance analysis: plotting the covariance matrix every T
    % timepoints after knocking out neurons cannot be grouped with other
    % neurons
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1.2.1  Test of algorithm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nskip             = false;
    if ~nskip
        nGroup        = 2;
        minSizeGroup  = 3;

        figure;
        for nPlot           = 1:numPlot
            subplot(m, m, nPlot)
            slicedDFF       = dff(:,timePoints(nPlot)+1:timePoints(nPlot+1));
            groupInd        = clusterWithKickOut(slicedDFF, nGroup, minSizeGroup);
            % regroupInd      = reGroupIndex(groupInd, nGroup);
            if sum(groupInd)>0
                % imagesc(corr(slicedDFF(regroupInd,:)'),[-1 1]);
                slicedDFF   = slicedDFF(groupInd>0, :);
                distNeurons     = pdist(slicedDFF, 'correlation');
                linkNeurons     = linkage(slicedDFF,'single','correlation');
                leafOrder       = optimalleaforder(linkNeurons, distNeurons);
                imagesc(corr(slicedDFF(leafOrder,:)'),[-1 1]);
                xlabel('Neuronal index')
                ylabel('Neuronal index')
                title(['Time from ' num2str((nPlot-1)*timeStep/1000) 'k to ' num2str(nPlot*timeStep/1000) 'k']);
                box off;
            end
        end
        setPrint(m*8, m*6, [plotDir, 'KnockingOut_UngroupedNeurons_', fileName], 'pdf');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1.2.2  Computing a minimal set of ungrouping neurons
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nskip                   = true;
    if ~nskip
        nGroup              = 2;
        minSizeGroup        = 3;
        numNeurons          = size(dff, 1);
        indGroup            = zeros(numNeurons, numPlot);

        for nPlot           = 1:numPlot
            slicedDFF       = dff(:,timePoints(nPlot)+1:timePoints(nPlot+1));
            groupInd        = clusterWithKickOut(slicedDFF, nGroup, minSizeGroup);
            indGroup(:,nPlot) = groupInd>0;
        end

        % 0: ungroup neurons           -- gray
        % 1: neurons in group          -- blue
        % 2: new neurons in group      -- red
        % 3: leaving-out neurons       -- purple

        cColor                   = [150 150 150; 23 100 171; 187 20 25; 23 171 25]/256;

        ind4Group                = zeros(numNeurons, numPlot);
        indUngroup               = sum(indGroup, 2) == 0;
        ind4Group(~indUngroup,:) = 1;
        for nPlot                = 1:numPlot - 1
            ind4Group(indGroup(:,nPlot)==0 & indGroup(:,nPlot+1)==1,nPlot+1) = 2;
            ind4Group(indGroup(:,nPlot+1)==0 & indGroup(:,nPlot)==1,nPlot+1) = 3;
        end

        % 2D plot    
    %     figure;
    %     for nPlot           = 1:numPlot
    %         subplot(m, m, nPlot)
    %         track_x         = squeeze(mean(tracks(:, timePoints(nPlot)+1:timePoints(nPlot+1),1), 2));
    %         track_y         = squeeze(mean(tracks(:, timePoints(nPlot)+1:timePoints(nPlot+1),2), 2));
    %         track_z         = squeeze(mean(tracks(:, timePoints(nPlot)+1:timePoints(nPlot+1),3), 2));
    %         var_x           = squeeze(var(tracks(:, timePoints(nPlot)+1:timePoints(nPlot+1),1), [], 2));
    %         var_y           = squeeze(var(tracks(:, timePoints(nPlot)+1:timePoints(nPlot+1),2), [], 2));
    %         var_z           = squeeze(var(tracks(:, timePoints(nPlot)+1:timePoints(nPlot+1),3), [], 2));
    %         radius          = sqrt(var_x + var_y + var_z);
    %         scatter3(track_x, track_y, track_z, radius, cColor(ind4Group(:,nPlot)+1,:));
    %         title(['Time from ' num2str(nPlot-1) 'k to ' num2str(nPlot) 'k']);
    %     end

        % 3D plot
        figure;
        for nPlot           = 1:numPlot
            subplot(m, m, nPlot)
            track_x         = squeeze(mean(tracks(:, timePoints(nPlot)+1:timePoints(nPlot+1),1), 2));
            track_y         = squeeze(mean(tracks(:, timePoints(nPlot)+1:timePoints(nPlot+1),2), 2));
            track_z         = squeeze(mean(tracks(:, timePoints(nPlot)+1:timePoints(nPlot+1),3), 2));
            var_x           = squeeze(var(tracks(:, timePoints(nPlot)+1:timePoints(nPlot+1),1), [], 2));
            var_y           = squeeze(var(tracks(:, timePoints(nPlot)+1:timePoints(nPlot+1),2), [], 2));
            var_z           = squeeze(var(tracks(:, timePoints(nPlot)+1:timePoints(nPlot+1),3), [], 2));
            radius          = sqrt(var_x + var_y + var_z);
            scatter3sph(track_x, track_y, track_z, 'siz', radius*5, 'col', cColor(ind4Group(:,nPlot)+1,:));
            title(['Time from ' num2str((nPlot-1)*timeStep/1000) 'k to ' num2str(nPlot*timeStep/1000) 'k']);
        end

        savefig([plotDir, 'Correlation_Clustering_', fileName]);
        save([tempDatDir, fileName, '.mat'], 'indUngroup', '-append');
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1.2.3  Using the last time point to estimate the group
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nskip             = true;
    if ~nskip
        nGroup        = 2;
        minSizeGroup  = 5;
        slicedDFF     = dff(:,timePoints(numPlot-1)+1:timePoints(numPlot));
        groupInd      = clusterWithKickOut(slicedDFF, nGroup, minSizeGroup);
        slicedDFF     = slicedDFF(groupInd>0, :);
        distNeurons   = pdist(slicedDFF, 'correlation');
        linkNeurons   = linkage(slicedDFF,'single','correlation');
        leafOrder     = optimalleaforder(linkNeurons, distNeurons);

        figure;
        for nPlot         = 1:numPlot
            subplot(m, m, nPlot)
            slicedDFF     = dff(:,timePoints(nPlot)+1:timePoints(nPlot+1));
            slicedDFF     = slicedDFF(groupInd>0, :);
            imagesc(corr(slicedDFF(leafOrder,:)'),[-1 1]);
            xlabel('Neuronal index')
            ylabel('Neuronal index')
            xlim([0 sum(groupInd>0)])
            set(gca,'Xtick',0:20:sum(groupInd>0))
            set(gca,'Ytick',0:20:sum(groupInd>0))
            title(['Time from ' num2str(floor(timePoints(nPlot)/1000)) 'k to ' num2str(floor(timePoints(nPlot+1)/1000)) 'k']);
            box off;
        end
        setPrint(m*8, m*6, [plotDir, 'Correlation_OderbyLastTimeSlice_', fileName], 'tif');
    end
    
end
