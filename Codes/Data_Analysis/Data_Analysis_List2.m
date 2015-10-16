%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  Covariance analysis: factor analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

function Data_Analysis_List2(nFile)
    % load data
    addpath('../Func');
    setDir;
    
    fileDirName       = fileDirNames{nFile}; %#ok<USENS,NASGU>
    fileName          = fileNames{nFile}; %#ok<USENS>
    
    % Data names:
    % dffWithNoDuplicates
    % indWithNoDuplicates
    % tracksWithNoDuplicates
    % indUngroup (after Data analysis #1)
    load([tempDatDir, fileName, '.mat']);
    %
    %
    % clear dffWithNoDuplicates indWithNoDuplicates indUngroup;
    maxNumFactor  = 20;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1.1  Covariance analysis: Factor analysis -- number of optimal factors
    % based on AIC or BIC
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1.1.1  Covariance analysis: Factor analysis -- number of optimal factors
    % based on AIC or BIC -- All neurons
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % timeStep      = 10000;
    % [numUnit, numT]  = size(dff);
    % numPlot       = ceil(numT/timeStep);
    % timePoints    = [0 (1:(numPlot-1))*timeStep numT];
    % 
    % matAIC        = zeros(maxNumFactor, numPlot);
    % matBIC        = zeros(maxNumFactor, numPlot);
    % matLL         = zeros(maxNumFactor, numPlot);
    % 
    % nskip             = true;
    % if ~nskip
    %     dff           = dffWithNoDuplicates;
    %     for nPlot        = 1:numPlot
    %         slicedDFF    = dff(:,timePoints(nPlot)+1:timePoints(nPlot+1));
    %         for nFactor  = 1:maxNumFactor
    %             [~, ~, ~, stats, ~]   = factoran(slicedDFF', nFactor, 'rotate', 'none');
    %             p_star                = numUnit * (nFactor+1) - (nFactor-1)*nFactor/2;
    %             n_star                = timePoints(nPlot+1) - timePoints(nPlot);
    %             matLL(nFactor,nPlot)  = stats.loglike;
    %             matAIC(nFactor,nPlot) = -2*stats.loglike + 2*p_star;
    %             matBIC(nFactor,nPlot) = -2*stats.loglike + log(n_star)*p_star;
    %         end
    %     end
    %     figure;
    %     subplot(1, 2, 1)
    %     imagesc(matAIC);
    %     title('AIC')
    %     subplot(1, 2, 2)
    %     imagesc(matBIC);
    %     title('BIC')
    %     suptitle('All neurons')
    % end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1.1.2  Covariance analysis: Factor analysis -- number of optimal factors
    % based on AIC or BIC -- grouped neurons
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % nskip             = true;
    % if ~nskip
    %     dff           = dffWithNoDuplicates(~indUngroup, :);
    %     for nPlot        = 1:numPlot
    %         slicedDFF    = dff(:,timePoints(nPlot)+1:timePoints(nPlot+1));
    %         for nFactor  = 1:maxNumFactor
    %             [~, ~, ~, stats, ~]   = factoran(slicedDFF', nFactor);
    %             p_star                = numUnit * (nFactor+1) - (nFactor-1)*nFactor/2;
    %             n_star                = timePoints(nPlot+1) - timePoints(nPlot);
    %             matLL(nFactor,nPlot)  = stats.loglike;
    %             matAIC(nFactor,nPlot) = -2*stats.loglike + 2*p_star;
    %             matBIC(nFactor,nPlot) = -2*stats.loglike + log(n_star)*p_star;
    %         end
    %     end
    %     figure;
    %     subplot(1, 2, 1)
    %     imagesc(matAIC);
    %     title('AIC')
    %     subplot(1, 2, 2)
    %     imagesc(matBIC);
    %     title('BIC')
    %     suptitle('Grouped neurons')
    % end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1.1.1  Covariance analysis: Factor analysis -- number of optimal factors
    % based on EV -- All neurons
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % timeStep      = 1200;
    % dff           = dffWithNoDuplicates;
    % [numUnit, numT]  = size(dff);
    % numPlot       = ceil(numT/timeStep);
    % timePoints    = [0 (1:(numPlot-1))*timeStep numT];
    numFold       = 10;
    numPlot       = length(timePoints)-1;
    matEV         = nan(numPlot, maxNumFactor, numFold);
    numUnit       = size(dff,1); %#ok<NODEF>
    matEVSingleUnit = nan(numPlot, maxNumFactor, numUnit, numFold);



    nskip             = false;
    if ~nskip    
    %     figure;

        for nPlot        = numPlot:-1:1
            slicedDFF    = dff(:,timePoints(nPlot)+1:timePoints(nPlot+1));
            slicedDFF    = bsxfun(@minus, slicedDFF, mean(slicedDFF,2));
            slicedDFF    = bsxfun(@rdivide, slicedDFF, std(slicedDFF,[],2));
            numSample    = size(slicedDFF, 2);
            numTest      = ceil(numSample * 0.1);
            for nFactor         = 1:maxNumFactor
                for nFold       = 1:numFold
                    randSeqTest                  = randperm(numSample);
                    nFoldDFFTest                 = slicedDFF(:,randSeqTest(1:numTest))';
                    nFoldDFFTrain                = slicedDFF(:,randSeqTest(numTest+1:numSample))';
                    [lambda,psi]                 = factoran(nFoldDFFTrain, nFactor, 'scores','regression', 'rotate', 'none');
                    matEV(nPlot, nFactor, nFold) = LONOFA(nFoldDFFTest, lambda, psi);
                    matEVSingleUnit(nPlot, nFactor, :, nFold) = LONOFASingleUnitEV(nFoldDFFTest, lambda, psi);
                    display([nPlot, nFactor, nFold])
                end
            end
    %         subplot(1,2,1);
    %         imagesc(timePoints, 1:maxNumFactor, mean(matEV,3)');
    %         subplot(1,2,2);
    %         plot(1:maxNumFactor, mean(matEV(nPlot,:,:),3)');        
    %         pause(0.1);
        end
    end

    % figure;
    % plot(1:maxNumFactor, mean(matEV,3))
    % errorbar(repmat(timePoints(1:end-1),maxNumFactor,1), mean(matEV,3)', std(matEV,[],3)');

    save([tempDatDir, fileName, '_FAEV.mat'], 'matEV');

    max_ev_units = squeeze(max(squeeze(mean(matEVSingleUnit,4)),[],2));
    figure ;
    for ii = 1:3
        subplot(1,3,ii);
        set(0,'DefaultAxesColorOrder',distinguishable_colors(20))
        plot(timePoints(1:end-1)/4/3600,max_ev_units(:,idx==ii)),ylim([0 1]) %#ok<COLND>
    end
    save([tempDatDir, fileName, '_FAEV.mat'], 'matEVSingleUnit','-append');
    
end