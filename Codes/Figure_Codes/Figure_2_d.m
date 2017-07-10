% Fit for Figure 2c for each fish
% 1. number of communities
% 2a. Fraction of non-factored neurons
% 2b. percentage of total active neurons
% 3a. radius of communities
% 3b. size of communities


function Figure_2_d(nFile)
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>  
    load([tempDatDir, 'FALONO_', fileName, '.mat'], 'EVLONO');
    load([tempDatDir, fileName, '.mat'],'activeNeuronMat', 'side', 'mnx');
    if ~exist([tempDatDir, 'EvoLoading_' fileName, '_v2.mat'], 'file'); return; end
    load([tempDatDir, 'EvoLoading_' fileName, '_v2.mat'], 'networkMat')
        
    subplot(1, 4, 1)
    Figure_2_d_1(EVLONO)    
    
    subplot(1, 4, 2)
    Figure_2_d_2(activeNeuronMat, networkMat)
end

%% 1. number of communities
function Figure_2_d_1(EVLONO)    
    numTime           = size(EVLONO, 1);
    LONOM             = zeros(numTime, 1);
    EVLONOMat         = squeeze(mean(EVLONO, 3));
    for nTime         = 1:numTime    
        if sum(~isnan(EVLONOMat(nTime, :))) > 0
            maxEV        = nanmax(EVLONOMat(nTime, :));
            if maxEV>0
                LONOM(nTime) = find(EVLONOMat(nTime, :)>0.9*maxEV, 1, 'first');
            end
        end        
    end
    timePoints        = (1:numTime)';
    hold on    
    plot(timePoints/60, LONOM,'o', 'color', [0.7 0.7 0.7])    
    % fit of numFactor curve
    fitResult    = fit(timePoints/60, LONOM-2, 'gauss1');
    b            = fitResult.b1;
    a            = fitResult.a1;
    c            = fitResult.c1;
    cr           = c;
    fitResult    = lsqcurvefit(@(p, x) doubleSizedGauss(p, x), [a, b, c, cr], timePoints/60, LONOM);    
    opt1Dim      = doubleSizedGauss(fitResult,timePoints/60);    
    plot(timePoints/60, opt1Dim,'k-', 'linewid', 2)
end

%% 2. 
% a. Fraction of non-factored neurons 
% b. percentage of total active neurons
function Figure_2_d_2(activeNeuronMat, networkMat)    
    
    numActNeuron      = sum(activeNeuronMat, 1);
    fracActNeuron     = mean(activeNeuronMat, 1);    
    numTime           = length(numActNeuron);
    numFactor         = zeros(numTime, 1);
    factorNeuronMat   = false(size(activeNeuronMat));
    % computer number of factored neuron
    for nTime         = 1:numTime
        factorSet     = networkMat{nTime, 1};
        for nFactor   = 2:length(factorSet)
            if length(factorSet(nFactor).neuronIndex) > 1
                factorNeuronMat(factorSet(nFactor).neuronIndex, nTime) = true;
                numFactor(nTime)       = numFactor(nTime)+1;
            end
        end
    end
    numFactorNeuron   = sum(factorNeuronMat, 1);
    fracNeuron        = numFactorNeuron./numActNeuron;
    
    hold on
    plot((1:numTime)/60, fracActNeuron, 'ok')
    plot((1:numTime)/60, 1 - fracNeuron, 'or')
    
    % fit for fracActNeuron
    % remove drop points
    removeIndex              = find(fracActNeuron/max(fracActNeuron) < 0.6);
    factorTime               = find(numFactor>= max(numFactor)-1, 1, 'last');
    removeIndex(removeIndex<factorTime) = [];
    if isempty(removeIndex)
        removeStart          = numTime;
    else
        removeStart          = removeIndex(1);
    end
    % for nfile = 12: removeStart  = 180; fix_fit_max = 0.95;
    % for nfile = 13; fix_fit_max = 0.95;
    % for nfile = 17: removeStart  = 180; fix_fit_max = 0.95;
    % for nfile = 18; removeStart  = 120;
    % for nfile = 19; removeStart  = 150; fix_fit_max = 0.95;
%     if nFile >=17; removeStart  = 180; end
    init_params              = [0, max(fracActNeuron), find(fracActNeuron/max(fracActNeuron)>0.5, 1, 'first')/60, 1];
    [fitParams, fitResult]   = sigm_fit((1:removeStart)/60, fracActNeuron(1:removeStart), [0, nan, nan, nan], init_params, false);
    upK                      = fitResult.paramCI(end, 1);
    lowK                     = fitResult.paramCI(end, 2);
    ypred                    = fitResult.ypred;
    ypredlowerCI             = fitResult.ypredlowerCI;
    ypredupperCI             = fitResult.ypredupperCI;
    plot((1:removeStart)/60, ypred, '-', 'linewid', 2.0, 'Color', 'k');
    plot((1:removeStart)/60, ypredlowerCI, '-', 'linewid', 0.5, 'Color', 'k');
    plot((1:removeStart)/60, ypredupperCI, '-', 'linewid', 0.5, 'Color', 'k');
    ylim([0 1])
    xlim([0 numTime/60])    
        
    % remove drop points
    removeIndex              = find(fracNeuron < 0.6);
    factorTime               = find(numFactor>= max(numFactor)-1, 1, 'last');
    removeIndex(removeIndex<factorTime) = [];
    if isempty(removeIndex)
        removeStart          = numTime;
    else
        removeStart          = removeIndex(1);
    end
    % for nfile = 12: removeStart  = 180; fix_fit_max = 0.95;
    % for nfile = 13; fix_fit_max = 0.95;
    % for nfile = 17: removeStart  = 180; fix_fit_max = 0.95;
    % for nfile = 18; removeStart  = 120;
    % for nfile = 19; removeStart  = 150; fix_fit_max = 0.95;
%     if nFile >=17; removeStart  = 180; end
    init_params              = [0, max(fracNeuron), find(fracNeuron>0.5, 1, 'first')/60, 1];
    [fitParams, fitResult]   = sigm_fit((1:removeStart)/60, fracNeuron(1:removeStart), [0, 0.95, nan, nan], init_params, false);
    upK                      = fitResult.paramCI(end, 1);
    lowK                     = fitResult.paramCI(end, 2);
    ypred                    = fitResult.ypred;
    ypredlowerCI             = fitResult.ypredlowerCI;
    ypredupperCI             = fitResult.ypredupperCI;
    plot((1:removeStart)/60, 1 - ypred, '-', 'linewid', 2.0, 'Color', 'r');
    plot((1:removeStart)/60, 1 - ypredlowerCI, '-', 'linewid', 0.5, 'Color', 'r');
    plot((1:removeStart)/60, 1- ypredupperCI, '-', 'linewid', 0.5, 'Color', 'r');
    ylim([0 1])
    xlim([0 numTime/60])
    
    
end

%% 3a. radius of communities
function Figure_2_d_3a(EVLONO)    
    numTime           = size(EVLONO, 1);
    LONOM             = zeros(numTime, 1);
    EVLONOMat         = squeeze(mean(EVLONO, 3));
    for nTime         = 1:numTime    
        if sum(~isnan(EVLONOMat(nTime, :))) > 0
            maxEV        = nanmax(EVLONOMat(nTime, :));
            if maxEV>0
                LONOM(nTime) = find(EVLONOMat(nTime, :)>0.9*maxEV, 1, 'first');
            end
        end        
    end
    timePoints        = (1:numTime)';
    hold on    
    plot(timePoints/60, LONOM,'o', 'color', [0.7 0.7 0.7])    
    % fit of numFactor curve
    fitResult    = fit(timePoints/60, LONOM-2, 'gauss1');
    b            = fitResult.b1;
    a            = fitResult.a1;
    c            = fitResult.c1;
    cr           = c;
    fitResult    = lsqcurvefit(@(p, x) doubleSizedGauss(p, x), [a, b, c, cr], timePoints/60, LONOM);    
    opt1Dim      = doubleSizedGauss(fitResult,timePoints/60);    
    plot(timePoints/60, opt1Dim,'k-', 'linewid', 2)
end

%% 3b. size of communities
function Figure_2_d_3b(EVLONO)    
    numTime           = size(EVLONO, 1);
    LONOM             = zeros(numTime, 1);
    EVLONOMat         = squeeze(mean(EVLONO, 3));
    for nTime         = 1:numTime    
        if sum(~isnan(EVLONOMat(nTime, :))) > 0
            maxEV        = nanmax(EVLONOMat(nTime, :));
            if maxEV>0
                LONOM(nTime) = find(EVLONOMat(nTime, :)>0.9*maxEV, 1, 'first');
            end
        end        
    end
    timePoints        = (1:numTime)';
    hold on    
    plot(timePoints/60, LONOM,'o', 'color', [0.7 0.7 0.7])    
    % fit of numFactor curve
    fitResult    = fit(timePoints/60, LONOM-2, 'gauss1');
    b            = fitResult.b1;
    a            = fitResult.a1;
    c            = fitResult.c1;
    cr           = c;
    fitResult    = lsqcurvefit(@(p, x) doubleSizedGauss(p, x), [a, b, c, cr], timePoints/60, LONOM);    
    opt1Dim      = doubleSizedGauss(fitResult,timePoints/60);    
    plot(timePoints/60, opt1Dim,'k-', 'linewid', 2)
end


%% related functions
function y = doubleSizedGauss(p, x)
    p      =  num2cell(p);
    [a, b, c, cr] = deal(p{:});
    d = 0;
    dr = 2;    
    y = x;
    y(x<=b) = a*exp(-((x(x<=b)-b)/c).^2) + d;
    y(x>b)  = (a+d-dr)*exp(-((x(x>b)-b)/cr).^2) + dr;
end