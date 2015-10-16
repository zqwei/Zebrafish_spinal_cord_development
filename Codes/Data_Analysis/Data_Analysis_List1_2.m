%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.  Covariance analysis: clustering analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

function Data_Analysis_List1_2(nFile)
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
    % Check results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %
    timeStep      = 1200;
    numT          = size(dff, 2);
    numPlot       = ceil(numT/timeStep);
    timePoints    = [0 (1:(numPlot-1))*timeStep numT];
    m             = ceil(sqrt(numPlot + 1));

    % nskip             = false;
    % if ~nskip  
    %     figure;
    %     for nPlot           = 1:numPlot
    %         subplot(m, m, nPlot)
    %         slicedDFF       = dff(:,timePoints(nPlot)+1:timePoints(nPlot+1));        
    %         imagesc(corr(slicedDFF'),[-1 1]);
    %         xlabel('Neuronal index')
    %         ylabel('Neuronal index')
    %         title(['Time from ' num2str((nPlot-1)*timeStep/1000) 'k to ' num2str(nPlot*timeStep/1000) 'k']);
    %         box off;
    %     end
    % end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1.2  Covariance analysis: plotting the covariance matrix every T
    % timepoints by removing the points with weak linear correlations with the
    % other points in the neuronal pool
    % KMO & Bartlett's test
    % KMO is for the sample adequacy
    % Bartlett's test checks if there is certain redundancy between the
    % variables that we can summarize with a few factors
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1.2.1 Bartlett's test of Sphericity, whether the correlation matrix is
    % identity matrix (Of note, barttest in matlab is for equal variance test)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    nskip             = true;
    ndimSlicedDFF     = zeros(numPlot-1,1);
    alpha             = 0.05;
    if ~nskip
    %     figure;
        for nPlot           = 1:numPlot-1
            slicedDFF       = dff(:,timePoints(nPlot)+1:timePoints(nPlot+1))';
            meanDFF         = mean(slicedDFF);
            stdDFF          = std(slicedDFF, 1);
            slicedDFF       = bsxfun(@rdivide, bsxfun(@minus, slicedDFF, meanDFF), stdDFF);
            Barspher(slicedDFF,alpha)
            ndimSlicedDFF(nPlot) = ndim;
        end    
    %     plot(timePoints(1:numPlot-1), ndimSlicedDFF,'o')   
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % All subsessions are factorable.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1.2.2 Correlation matrix -- number of eigenvalues > 1.0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nskip                   = false;
    logDetCorrSlicedDFF     = zeros(numPlot-1,1);
    logDetCorrAllDFF        = zeros(numPlot-1,1);
    % thres                   = 1.0;
    if ~nskip
        figure;
        for nPlot           = 1:numPlot-1
            allDFF          = dffWithNoDuplicates(:,timePoints(nPlot)+1:timePoints(nPlot+1))';
            slicedDFF       = dff(:,timePoints(nPlot)+1:timePoints(nPlot+1))';
            % ndimSlicedDFF(nPlot) = sum(eig(corr(slicedDFF))>thres);
            logDetCorrSlicedDFF(nPlot) = log(det(corr(slicedDFF)));
            logDetCorrAllDFF(nPlot)    = log(det(corr(allDFF)));
        end    

        timePoints          = timePoints(logDetCorrAllDFF>=logDetCorrAllDFF(end));
        logDetCorrSlicedDFF = logDetCorrSlicedDFF(logDetCorrAllDFF>=logDetCorrAllDFF(end));
        plot(timePoints/4/3600, logDetCorrSlicedDFF,'o')  %ndimSlicedDFF/size(slicedDFF,1)
        ylabel('log |corr(dF/F)|')
        xlabel('time (hours)')
        box off
        setPrint(8, 6, [plotDir 'LogDetCorrAllNeurons_' fileName], 'pdf')
    end

    save([tempDatDir, fileName, '.mat'],'timePoints','-append')


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1.2.3 Correlation matrix -- KMO -- decides whether sample number is
    % enough in the 0.90 as marvelous, in the 0.80's as meritorious, in the 
    % 0.70's as middling, in the 0.60's as mediocre, in the 0.50's as miserable
    % , and below 0.50 as unacceptable
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nskip             = true;
    unitKMO           = zeros(numPlot-1,size(dff,1));
    totKMO            = zeros(numPlot-1,1);
    if ~nskip
        for nPlot                 = 1:numPlot-1
            slicedDFF             = dff(:,timePoints(nPlot)+1:timePoints(nPlot+1))';
            [totKMO(nPlot), unitKMO(nPlot,:)] = kmo(slicedDFF);
        end
        figure;
        imagesc(timePoints(1:numPlot-1), 1:size(dff,1), unitKMO);
        caxis([0.5 1.0])
        figure;
        plot(timePoints(1:numPlot-1), totKMO,'o')
    end
    
end

