%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  Covariance analysis: factor analysis
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.1 Factor as a function of the time and dimensions of factors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

function Data_Analysis_List2_1(nFile)
    addpath('../Func');
    setDir;    
    fileDirName       = fileDirNames{nFile};
    fileName          = fileNames{nFile};
    
    load([tempDatDir, fileName, '.mat']);
    maxNumFactor  = 20;

    load([tempDatDir, fileName, '_FAEV.mat'], 'matEV', 'matEVSingleUnit');
    figure;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   plot EV as a function of time and Dim
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(2,3,[1 2]);
    plot_EV_Time_Dim(timePoints, maxNumFactor, matEV);
    axis xy;
    colorbar
    caxis([0 0.7])
    box off;
    title('Total EV')
    xlabel('Time (hour)')
    ylabel('Latent dimensions')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   plot optimal dimensionality for each time point
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subplot(2,3,3)
    opt1Dim = plot_optDim_expFit(timePoints, matEV, maxNumFactor);
    opt1Dim = round(opt1Dim);
    save([tempDatDir, fileName, '_FAEV.mat'], 'opt1Dim', '-append');
    opt2Dim = plot_optDim_expFitSmooth(timePoints, matEV, maxNumFactor);
    opt2Dim = round(opt2Dim); 
    if exist('opt2Dim','var'); opt1Dim = opt2Dim; end
    save([tempDatDir, fileName, '_FAEV.mat'], 'opt2Dim', '-append');
    box off;
    xlabel('Time (hour)')
    ylabel('Optimal Latent dimensions')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   plot EV for each single unit using optimal dimensions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    max_ev_units        = nan(size(matEVSingleUnit,3), length(timePoints)-1); %#ok<NODEF>
    for nTime           = 1:length(timePoints)-1
        max_ev_units(:,nTime) = squeeze(mean(matEVSingleUnit(nTime,opt1Dim(nTime),:,:),4));
    end
    subplot(2,3,4)
    hold on;
    plot(timePoints(1:end-1)/4/3600,max_ev_units(idx==1,:),'-','Color',[0.5 0.5 0.5],'linewid',0.5); %#ok<COLND>
    plot(timePoints(1:end-1)/4/3600,mean(max_ev_units(idx==1,:),1),'-k','linewid',1.0); %#ok<COLND>
    ylim([0 1]);
    xlim([0 timePoints(end)'/4/3600]) %#ok<COLND>
    title('Moving left neurons');
    xlabel('Time (hour)')
    ylabel('EV for each neuron')
    
    subplot(2,3,5)
    hold on;
    plot(timePoints(1:end-1)/4/3600,max_ev_units(idx==2,:),'-','Color',[0.5 0.5 0.5],'linewid',0.5); %#ok<COLND>
    plot(timePoints(1:end-1)/4/3600,mean(max_ev_units(idx==2,:),1),'-k','linewid',1.0); %#ok<COLND>
    ylim([0 1]);
    xlim([0 timePoints(end)'/4/3600]) %#ok<COLND>
    title('Moving right neurons');
    xlabel('Time (hour)')
    ylabel('EV for each neuron')
    
    subplot(2,3,6)
    hold on;
    plot(timePoints(1:end-1)/4/3600,max_ev_units(idx==3,:),'-','Color',[0.5 0.5 0.5],'linewid',0.5); %#ok<COLND>
    plot(timePoints(1:end-1)/4/3600,mean(max_ev_units(idx==3,:),1),'-k','linewid',1.0); %#ok<COLND>
    ylim([0 1]);
    xlim([0 timePoints(end)'/4/3600]) %#ok<COLND>
    title('Non-selective neurons');
    xlabel('Time (hour)')
    ylabel('EV for each neuron')
    
%     setPrint(8*3, 6*2, [plotDir, 'FACrossTime_', fileName], 'tif');
    
end

function plot_EV_Time_Dim(timePoints, maxNumFactor, matEV)
    imagesc(timePoints(1:end-1)/4/3600, 1:maxNumFactor, mean(matEV,3)');
    xlim([0 timePoints(end)'/4/3600])
    ylim([1 20])
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


function opt1Dim = plot_optDim_expFitSmooth(timePoints, matEV, maxNumFactor)
    opt1Dim      = plot_optDim_expFit(timePoints, matEV, maxNumFactor);
    plot(timePoints(1:end-1)'/4/3600,round(opt1Dim),'o','linewid',2);
    fitResult    = fit(timePoints(1:end-1)'/4/3600, opt1Dim - min(opt1Dim), 'gauss1');
    opt1Dim      = feval(fitResult,timePoints(1:end-1)'/4/3600) + min(opt1Dim);
    hold on;
    plot(timePoints(1:end-1)'/4/3600,opt1Dim,'-','linewid',2);
    hold off;
end

function opt1Dim = plot_optDim_anova1(timePoints, matEV)
    opt1Dim             = nan(length(timePoints)-1,1);% p-value = 0.01
    opt2Dim             = nan(length(timePoints)-1,1);% p-value = 0.05

    for nTime           = 1:length(timePoints)-1
        EVTime          = squeeze(mean(matEV(nTime,:,:),3));
        [~, maxDim]     = max(EVTime);
        p               = ones(maxDim,1);
        for nDim        = 1:maxDim-1
            p(nDim)     = anova1(squeeze(matEV(nTime,nDim:maxDim,:))',[],'off');
        end
        opt1Dim (nTime) = find(p>0.001, 1, 'first');
        opt2Dim (nTime) = find(p>0.05, 1, 'first');
    end
    plot(timePoints(1:end-1)'/4/3600,opt1Dim,'o','linewid',2);
    xlim([0 timePoints(end)'/4/3600])
end