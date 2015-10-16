%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  Covariance analysis: factor analysis
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.0 Figure summary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

function Data_Analysis_List2_0(nFile)
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>
    
    load([tempDatDir, fileName, '.mat']);
    maxNumFactor      = 20;
    maxNumFactor2     = 10; 
    minNumFactor      = 2;
    EVSingleTime      = load([tempDatDir, fileName, '_FAEV.mat'], 'matEV', 'matEVSingleUnit', 'opt2Dim', 'opt1Dim');
    EVCrossTime       = load([tempDatDir, fileName, '_FAEVCrossTime.mat'], 'matEV', 'matEVSingleUnit');
    
    
    
    figure;   
%     h = suptitle(fileName);
%     set(h,'Interpreter','none');
    subplot(4,3,[1 2]);
    plot_EV_Time_Dim(timePoints, maxNumFactor, EVSingleTime.matEV);

    subplot(4,3,3)
    plot_optDim_expFitSmooth(timePoints, EVSingleTime.matEV, maxNumFactor)
    
    subplot(4,3,[4 5]);
    plot_EV_Time_Time(timePoints, EVCrossTime.matEV)
    
    subplot(4,3,[7 8]);
    opt1Dim = EVSingleTime.opt2Dim;
    opt1Dim = min(opt1Dim,maxNumFactor2-minNumFactor+1);
    plot_EVopt1Dim_Time_Time(timePoints, EVCrossTime.matEV, opt1Dim, minNumFactor)
    
    max_ev_units        = nan(size(EVCrossTime.matEVSingleUnit,4), length(timePoints)-1);
    for nTime           = 1:length(timePoints)-1
        max_ev_units(:,nTime) = squeeze(EVCrossTime.matEVSingleUnit(nTime,nTime,opt1Dim(nTime)-minNumFactor+1,:));
    end    
    
    subplot(4,3,10)
    hold on;
    plot(timePoints(1:end-1)/4/3600,max_ev_units(idx==1,:),'-','Color',[0.5 0.5 0.5],'linewid',0.5); %#ok<COLND>
    plot(timePoints(1:end-1)/4/3600,mean(max_ev_units(idx==1,:),1),'-k','linewid',1.0); %#ok<COLND>
    ylim([0 1]);
    xlim([0 timePoints(end)'/4/3600]) %#ok<COLND>
    title('Moving left neurons');
    xlabel('Time (hour)')
    ylabel('EV for each neuron')
    
    subplot(4,3,11)
    hold on;
    plot(timePoints(1:end-1)/4/3600,max_ev_units(idx==2,:),'-','Color',[0.5 0.5 0.5],'linewid',0.5); %#ok<COLND>
    plot(timePoints(1:end-1)/4/3600,mean(max_ev_units(idx==2,:),1),'-k','linewid',1.0); %#ok<COLND>
    ylim([0 1]);
    xlim([0 timePoints(end)'/4/3600]) %#ok<COLND>
    title('Moving right neurons');
    xlabel('Time (hour)')
    ylabel('EV for each neuron')
    
    subplot(4,3,12)
    hold on;
    plot(timePoints(1:end-1)/4/3600,max_ev_units(idx==3,:),'-','Color',[0.5 0.5 0.5],'linewid',0.5); %#ok<COLND>
    plot(timePoints(1:end-1)/4/3600,mean(max_ev_units(idx==3,:),1),'-k','linewid',1.0); %#ok<COLND>
    ylim([0 1]);
    xlim([0 timePoints(end)'/4/3600]) %#ok<COLND>
    title('Non-selective neurons');
    xlabel('Time (hour)')
    ylabel('EV for each neuron')

    % example of group one neuron
    subplot(4,3,6)
    [nNeuron, nTime]  = find(max_ev_units == max(max(max_ev_units(idx==1,:))));
    refDFF            = dff(:,timePoints(nTime)+1:timePoints(nTime+1)); %#ok<NODEF>
    refDFF            = bsxfun(@minus, refDFF, mean(refDFF,2));
    refDFF            = bsxfun(@rdivide, refDFF, std(refDFF,[],2))';
    [lambda,psi]      = factoran(refDFF, EVSingleTime.opt2Dim(nTime), 'scores','regression', 'rotate', 'none');
    DFFEst            = LONOFASingleUnit (refDFF, lambda, psi, nNeuron);
    hold on;
    plot((timePoints(nTime)+1:timePoints(nTime+1))/4/3600,refDFF(:,nNeuron),'-','linewid',0.5);
    plot((timePoints(nTime)+1:timePoints(nTime+1))/4/3600,DFFEst,'-','linewid',0.5);
    hold off
    xlim([(timePoints(nTime)+1)/4/3600 timePoints(nTime+1)/4/3600])
    title('Selective neurons');
    xlabel('Time (hour)')
    ylabel('Normalized DF/F')
    
    subplot(4,3,9)
    [nNeuron, nTime]  = find(max_ev_units == max(max(max_ev_units(idx==3,:))));
    refDFF            = dff(:,timePoints(nTime)+1:timePoints(nTime+1));
    refDFF            = bsxfun(@minus, refDFF, mean(refDFF,2));
    refDFF            = bsxfun(@rdivide, refDFF, std(refDFF,[],2))';
    [lambda,psi]      = factoran(refDFF, EVSingleTime.opt2Dim(nTime), 'scores','regression', 'rotate', 'none');
    DFFEst            = LONOFASingleUnit (refDFF, lambda, psi, nNeuron);
    hold on;
    plot((timePoints(nTime)+1:timePoints(nTime+1))/4/3600,refDFF(:,nNeuron),'-','linewid',0.5);
    plot((timePoints(nTime)+1:timePoints(nTime+1))/4/3600,DFFEst,'-','linewid',0.5);
    hold off
    xlim([(timePoints(nTime)+1)/4/3600 timePoints(nTime+1)/4/3600])
    title('Non-Selective neurons');
    xlabel('Time (hour)')
    ylabel('Normalized DF/F')
    setFigureSize(8*3, 6*4);
end

function plot_EV_Time_Dim(timePoints, maxNumFactor, matEV)
    imagesc(timePoints(1:end-1)/4/3600, 1:maxNumFactor, mean(matEV,3)');
    xlim([0 timePoints(end)'/4/3600])
    ylim([1 maxNumFactor])
    axis xy;
    colorbar;
    caxis([0 0.7])
    box off;
    title('Total EV')
    xlabel('Time (hour)')
    ylabel('Latent dimensions')
end


function opt1Dim = plot_optDim_expFit(timePoints, matEV, maxNumFactor)
    opt1Dim            = nan(length(timePoints)-1,1);
    fitFun             = 'c-exp(-(x-a)/b)';    
    approchThres       = 0.05;
    for nTime          = 1:length(timePoints)-1
        EVTime         = squeeze(mean(matEV(nTime,:,:),3));
        startPoint     = [0, 2 , max(EVTime)+0.001];
        fitResult      = fit((1:maxNumFactor)',EVTime',fitFun,'Start', startPoint);
        opt1Dim(nTime) = fitResult.a - fitResult.b*log(approchThres);
    end
    % opt1Dim            = round(opt1Dim);
    plot(timePoints(1:end-1)'/4/3600,round(opt1Dim),'o','linewid',2);
    xlim([0 timePoints(end)'/4/3600])
end


function plot_optDim_expFitSmooth(timePoints, matEV, maxNumFactor)
    opt1Dim      = plot_optDim_expFit(timePoints, matEV, maxNumFactor);
    fitResult    = fit(timePoints(1:end-1)'/4/3600, opt1Dim - min(opt1Dim), 'gauss1');
    opt1Dim      = feval(fitResult,timePoints(1:end-1)'/4/3600) + min(opt1Dim);
    hold on;
    plot(timePoints(1:end-1)'/4/3600,opt1Dim,'-','linewid',2);
    hold off;
    box off;
    xlabel('Time (hour)')
    ylabel('Smoothed Latent dimensions')
end


function plot_EV_Time_Time(timePoints, matEV)
    imagesc(timePoints(1:end-1)/4/3600, timePoints(1:end-1)/4/3600, max(matEV,[],3)');
    xlim([0 timePoints(end)'/4/3600])
    ylim([0 timePoints(end)'/4/3600])
    box off
    title ('Max Total EV')
    xlabel('Time (hour)') % reference
    ylabel('Time (hour)') % test
    colorbar
    caxis([-0.3 0.7])
end


function plot_EVopt1Dim_Time_Time(timePoints, matEV, opt1Dim, minNumFactor)
    maxEV      = zeros(size(matEV,1), size(matEV,2));
    
    for  nPlot = 1:length(timePoints)-1
        maxEV(:,nPlot) = squeeze(matEV(nPlot, :, opt1Dim(nPlot)-minNumFactor+1));
    end
    
    imagesc(timePoints(1:end-1)/4/3600, timePoints(1:end-1)/4/3600, maxEV);
    xlim([0 timePoints(end)'/4/3600])
    ylim([0 timePoints(end)'/4/3600])
    box off
    title ('Smoothed Total EV')
    xlabel('Time (hour)') % reference
    ylabel('Time (hour)') % test
    colorbar
    caxis([-0.3 0.7])
end


function DFFEst = LONOFASinglueUnit (nFoldDFFTest, lambda, psi, nUnit)


    LONODFF      = nFoldDFFTest(:, [1:nUnit-1 nUnit+1:end]);
    L            = lambda([1:nUnit-1 nUnit+1:end], :); % yDim -1 x xDim
    Ph           = psi([1:nUnit-1 nUnit+1:end]);
    
    Currlambda   = lambda(nUnit, :);
    
    estX         = L'/(L*L'+diag(Ph)) * LONODFF';
    % xDim x n
    % This method of estX is equivalent to F from following code
    % [lambda,psi, ~, ~, F] = factoran(data, m ,'scores','regression');
    
    DFFEst       = (Currlambda * estX)'; 
    

end