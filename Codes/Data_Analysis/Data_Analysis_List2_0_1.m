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

function Data_Analysis_List2_0_1(nFile)
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>
    
    load([tempDatDir, fileName, '.mat']);
    maxNumFactor      = 20;
    maxNumFactor2     = 10; 
    minNumFactor      = 2;
    load([tempDatDir, fileName, '_FAEV.mat'], 'opt2Dim');
    load([tempDatDir, fileName, '_FAEVCrossTime.mat'], 'matEV', 'matEVSingleUnit');
    
    opt1Dim          = opt2Dim-minNumFactor+1;
    opt1Dim          = min(opt1Dim, 8);
    
    numPlot          = length(timePoints)-1;
    mCol             = 8;
    mRow             = ceil((numPlot+1)/mCol);
    
    figure;
    h = suptitle(fileName);
    set(h,'Interpreter','none');
    for nPlot     = 1:numPlot
        subplot(mCol, mRow, nPlot)
        slicedDFF = dff(:,timePoints(nPlot)+1:timePoints(nPlot+1));
        imagesc(corr(slicedDFF'),[-1 1]);
        xlabel('Neuronal index')
        ylabel('Neuronal index')
        title(['Time from ' num2str(timePoints(nPlot)/3600/4) 'k to ' num2str(timePoints(nPlot+1)/3600/4) 'k']);
        box off;
    end
    % across time points
    subplot(mCol, mRow, numPlot + 1)
    slicedDFF = dff;
    imagesc(corr(slicedDFF'),[-1 1]);
    xlabel('Neuronal index')
    ylabel('Neuronal index')
    title(['Time from ' num2str(timePoints(1)/3600/4) 'k to ' num2str(timePoints(end)/3600/4) 'k']);
    box off;
    
    
    figure;   
    h = suptitle(fileName);
    set(h,'Interpreter','none');

    max_ev_units        = nan(size(matEVSingleUnit,4), length(timePoints)-1); %#ok<NODEF>
    for nTime           = 1:length(timePoints)-1
        max_ev_units(:,nTime) = squeeze(matEVSingleUnit(nTime,nTime,opt1Dim(nTime),:));
    end   
    
    numPlot          = size(max_ev_units,1);
    mCol             = 8;
    mRow             = ceil(numPlot/mCol);
    
    
    for nUnit        = 1:numPlot
        subplot(mRow, mCol, nUnit)
        if idx(nUnit) == 1
            mColor   = [ 0    0.4470    0.7410];
        elseif idx(nUnit) == 2
            mColor   = [0.6350    0.0780    0.1840];
        elseif idx(nUnit) == 3
            mColor   = [0 0 0];
        end
        plot(timePoints(1:end-1)/4/3600, max_ev_units(nUnit,:),'o','color',[0.5 0.5 0.5]); %#ok<COLND>
        
%         fitFun             = 'd+c/(1+exp(-(x-a)/b))';
%         startPoint         = [0, timePoints(end)/4/3600/3, max(max_ev_units(nUnit,:))-min(max_ev_units(nUnit,:)), min(max_ev_units(nUnit,:))]; %#ok<COLND>
%         fitResult          = fit((timePoints(1:end-1)/4/3600)', max_ev_units(nUnit,:)',fitFun,'Start', startPoint);
        [~, fitResult]     = sigm_fit((timePoints(1:end-1)/4/3600), max_ev_units(nUnit,:));
        hold on;
%         h                  = plot(fitResult);
%         h.LineStyle        = '-';
%         h.LineWidth        = 2.0;
%         h.Color            = mColor;
%         legend('off')
        plot((timePoints(1:end-1)/4/3600), fitResult.ypred, '-', 'linewid', 2.0, 'Color', mColor); %#ok<COLND>
        plot((timePoints(1:end-1)/4/3600), fitResult.ypredlowerCI, '-', 'linewid', 0.5, 'Color', mColor); %#ok<COLND>
        plot((timePoints(1:end-1)/4/3600), fitResult.ypredupperCI, '-', 'linewid', 0.5, 'Color', mColor); %#ok<COLND>
        hold off;
        
        RSquare             = 1 - mean((max_ev_units(nUnit,:)' - fitResult.ypred).^2)./var(max_ev_units(nUnit,:));        
        
        ylim([0 1])
        xlim([timePoints(1)/4/3600 timePoints(end)/4/3600]) %#ok<COLND>
%         title(num2str(nUnit,'%02d'))

        title(['R^2 = ', num2str(RSquare,'%.3f')])
        xlabel('Time (hour)')
        ylabel('EV')
        box off
    end
    
    setPrint(mCol*6, mRow*4.5, [plotDir, 'FAEVSingleUnit_', fileName], 'pdf');

end