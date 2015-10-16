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

function Data_Analysis_List2_1_0(nFile)
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>
    
    load([tempDatDir, fileName, '.mat'], 'timePoints');
    EVSingleTime      = load([tempDatDir, fileName, '_PSDPeakTimeFAEV.mat'], 'matEV');
    maxNumFactor      = 10;
    load([tempDatDir, fileName, '_numFactorNoLONOActiveNeuronsSimplified.mat'], 'numFactors', 'nActiveUnit');  
    
    figure;   
    subplot(2, 4,1:2);
    plot_EV_Time_Dim(timePoints, maxNumFactor, EVSingleTime.matEV);

    
    subplot(2, 4, 3);
    plot(timePoints(1:end-1)'/4/3600,nActiveUnit,'o','linewid',2);
    box off;
    xlim([0 timePoints(end)'/4/3600]);
    xlabel('Time (hour)')
    ylabel('# Active Neurons')

    
    fieldNames = {  'kgM', 'SRMRM', 'CFIM'};
    subplot(2, 4, 4+1)
    fitY         = zeros(length(timePoints)-1, 4);
    fitY(:, 1)   = plot_optDim_expFit(timePoints, EVSingleTime.matEV, maxNumFactor);
    title('LONO')
    xlabel('Time (hour)')
    ylabel('Latent dimensions')
    ylim([0 8])

%     hold on;
    for nPlot = 1:length(fieldNames)
        subplot(2, 4, 4+1+nPlot)
        numFactorCrossTime = nan(length(timePoints) - 1, 1);
        numIdVec = ~cellfun(@isempty, {numFactors.(fieldNames{nPlot})});
        numFactorCrossTime(numIdVec) = [numFactors.(fieldNames{nPlot})]';
        fitY(:, nPlot + 1) = numFactorCrossTime;
        plot(timePoints(1:end-1)'/4/3600,numFactorCrossTime,'o');
        title(fieldNames{nPlot})
        xlabel('Time (hour)')
        ylabel('Latent dimensions')
        ylim([0 8])
    end
    
%     opt2Dim   = round(mean(fitY(:, [1 3 4]), 2));
%     plot(timePoints(1:end-1)'/4/3600,opt2Dim,'xk');
%     ylim([0 maxNumFactor])
%     hold off
%     legend({'LONO', 'KG', 'SRMR', 'CFI', 'Fit'})    
%     legend('location','northeast')
%     legend('boxoff')
%     box off;
%     xlabel('Time (hour)')
%     ylabel('Latent dimensions')
    
    setFigureSize(8*4, 6*2);
    
%     save([tempDatDir, fileName, '_PSDPeakTimeFAEV.mat'], 'opt2Dim', '-append')

%     setPrint(8*3, 6*4, [plotDir, 'FACrossTime_', fileName], 'pdf');

end

function plot_EV_Time_Dim(timePoints, maxNumFactor, matEV)
    imagesc(timePoints(1:end-1)/4/3600, 1:maxNumFactor, nanmean(matEV,3)');
    xlim([0 timePoints(end)'/4/3600])
    ylim([1 maxNumFactor])
    axis xy;
    colorbar;
    caxis([0 0.8])
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
        EVTime         = squeeze(nanmean(matEV(nTime,:,:),3));
        if isnan(nanmax(EVTime))
            opt1Dim(nTime) = 0;
        else
            startPoint     = [0, 2 , nanmax(EVTime)+0.001];
            vecX           = (1:maxNumFactor)';
            vecY           = EVTime';
            vecX           = vecX(vecY>0);
            vecY           = vecY(vecY>0);
            if length(vecY) > 4
                fitResult      = fit(vecX, vecY,fitFun,'Start', startPoint);
                opt1Dim(nTime) = max(fitResult.a - fitResult.b*log(approchThres), 0);
            else
                [~, opt1Dim(nTime)] = max(vecX);
            end
        end
    end
    % opt1Dim            = round(opt1Dim);
    plot(timePoints(1:end-1)'/4/3600,round(opt1Dim),'o');
    xlim([0 timePoints(end)'/4/3600])
end


function plot_optDim_expFitSmooth(timePoints, matEV, maxNumFactor)
    opt1Dim      = plot_optDim_expFit(timePoints, matEV, maxNumFactor);
%     fitResult    = fit(timePoints(1:end-1)'/4/3600, opt1Dim - min(opt1Dim), 'gauss1');
%     opt1Dim      = feval(fitResult,timePoints(1:end-1)'/4/3600) + min(opt1Dim);
%     hold on;
%     plot(timePoints(1:end-1)'/4/3600,opt1Dim,'-','linewid',2);
%     hold off;
    box off;
    xlabel('Time (hour)')
    ylabel('Latent dimensions')
end


function plot_noLONODim_expFitSmooth(timePoints, numFactors, ylabelName)
%     fitResult    = fit(timePoints(1:end-1)'/4/3600, numFactors - min(numFactors), 'gauss1');
%     opt1Dim      = feval(fitResult,timePoints(1:end-1)'/4/3600) + min(numFactors);
    plot(timePoints(1:end-1)'/4/3600,numFactors,'o','linewid',2);
%     hold on;
%     plot(timePoints(1:end-1)'/4/3600,opt1Dim,'-','linewid',2);
%     hold off;
    box off;
    xlim([0 timePoints(end)'/4/3600]);
    xlabel('Time (hour)')
    ylabel(ylabelName)
end


