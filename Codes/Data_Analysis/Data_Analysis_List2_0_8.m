%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  Covariance analysis: factor analysis
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.0.3 Figure summary of number of factors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

function Data_Analysis_List2_0_8(nFile)
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>
    
    load([tempDatDir, fileName, '.mat'], 'timePoints');
                    
    load([tempDatDir, fileName, '_numFactorNoLONOActiveNeuronsSimplified.mat'], 'numFactors', 'nActiveUnit');                
        
    figure;       
    subplot(1, 4,1);
    plot(timePoints(1:end-1)'/4/3600,nActiveUnit,'o','linewid',2);
    box off;
    xlim([0 timePoints(end)'/4/3600]);
    xlabel('Time (hour)')
    ylabel('# Active Neurons')
        
    fieldNames = {  'kgM', 'SRMRM', 'CFIM'};
    ylabelNames = {  'Kaiser-Guttman Dim', 'SRMR Dim', 'CFI Dim'};

    for nPlot = 1:length(fieldNames)
        subplot(1, 4, 1+nPlot);
        numFactorCrossTime = nan(length(timePoints) - 1, 1);
        numIdVec = ~cellfun(@isempty, {numFactors.(fieldNames{nPlot})});
        numFactorCrossTime(numIdVec) = [numFactors.(fieldNames{nPlot})]';
        plot_noLONODim_expFitSmooth(timePoints, numFactorCrossTime, ylabelNames{nPlot});
    end

    
    setFigureSize(8*4, 6);
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


function plot_optDim_expFitSmooth(timePoints, matEV, maxNumFactor, ylabelName)
    opt1Dim      = plot_optDim_expFit(timePoints, matEV, maxNumFactor);
    fitResult    = fit(timePoints(1:end-1)'/4/3600, opt1Dim - min(opt1Dim), 'gauss1');
    opt1Dim      = feval(fitResult,timePoints(1:end-1)'/4/3600) + min(opt1Dim);
    hold on;
    plot(timePoints(1:end-1)'/4/3600,opt1Dim,'-','linewid',2);
    hold off;
    box off;
    xlabel('Time (hour)')
    ylabel(ylabelName)
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


