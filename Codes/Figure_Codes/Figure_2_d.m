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
    load([tempDatDir, 'EvoLoading_' fileName, '_v2.mat'], 'networkMat', 'neuronXLoc', 'neuronYLoc', 'neuronZLoc')
    load([tempDatDir, 'EV_' fileName, '.mat'], 'halfActTime', 'validFitIndex')
    
    numTime           = length(networkMat); %#ok<*USENS>
    clusterList       = []; % time, size, radius    
    for nTime         = 1:numTime
        factorSet     = networkMat{nTime, 1}; 
        for nFactor   = 2:length(factorSet)-1
            neuronIndex = factorSet(nFactor).neuronIndex;
            if length(neuronIndex)>1
                cluster = [ nTime,...
                            length(neuronIndex)];
                xLoc    = neuronXLoc(neuronIndex, nTime);
                yLoc    = neuronYLoc(neuronIndex, nTime);
                zLoc    = neuronZLoc(neuronIndex, nTime);
%                 radius  = (xLoc - factorSet(nFactor).x).^2 + (yLoc - factorSet(nFactor).y).^2 + ...
%                           (zLoc - factorSet(nFactor).z).^2;
                radius  = (xLoc - factorSet(nFactor).x).^2;
                radius  = sqrt(mean(radius));
                cluster = [cluster, radius];
                clusterList = [clusterList; cluster];
            end
        end
    end
    
    totPlots = 6;
        
    subplot(1, totPlots, 1)
    Figure_2_d_1(EVLONO)    
    
    subplot(1, totPlots, 2)
    Figure_2_d_2a(activeNeuronMat, networkMat)
    
    subplot(1, totPlots, 3)
    Figure_2_d_2b(activeNeuronMat, networkMat)
    
    subplot(1, totPlots, 4)
    Figure_2_d_3a(clusterList, numTime)
    
    subplot(1, totPlots, 5)
    Figure_2_d_3b(clusterList, numTime)
    
    subplot(1, totPlots, 6)
    Figure_2_d_4(halfActTime, neuronXLoc(:, 1))
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
    fitResult    = fit(timePoints/60, LONOM-2, 'gauss1');
    b            = fitResult.b1;
    a            = fitResult.a1;
    c            = fitResult.c1;
    cr           = c;
    fitResult    = lsqcurvefit(@(p, x) doubleSizedGauss(p, x), [a, b, c, cr], timePoints/60, LONOM);    
    opt1Dim      = doubleSizedGauss(fitResult,timePoints/60);  
    
    plot(timePoints/60, opt1Dim)
%     plotAlignPerc(timePoints/60, opt1Dim, 98)
    ylabel('Num factor')
    xlabel('Time from peak (hour)')
    box off
end

%% 2. 
% a. Fraction of non-factored neurons 
% b. percentage of total active neurons
function Figure_2_d_2a(activeNeuronMat, networkMat)    
    numActNeuron      = sum(activeNeuronMat, 1);
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
    % fit for fracNeuron    
    factorTime               = find(numFactor>= max(numFactor)-1, 1, 'last');
    fitResult                = fig_sigm((1:numTime)/60, fracNeuron, factorTime, 0, 1);
    timePoints               = fitResult.x;
    ypred                    = fitResult.ypred;
%     plotAlignPerc(timePoints, ypred, 90); 
    plot(timePoints, ypred); 
    box off
    ylabel('Frac factored neuron')
    xlabel('Time from 90% (hour)')
end

function Figure_2_d_2b(activeNeuronMat, networkMat)    
    numActNeuron      = sum(activeNeuronMat, 1);
    fracActNeuron     = mean(activeNeuronMat, 1);    
    numTime           = length(numActNeuron);
    numFactor         = zeros(numTime, 1);   
    for nTime         = 1:numTime
        factorSet     = networkMat{nTime, 1};
        for nFactor   = 2:length(factorSet)
            if length(factorSet(nFactor).neuronIndex) > 1
                numFactor(nTime)       = numFactor(nTime)+1;
            end
        end
    end
    hold on
    % fit for fracActNeuron    
    factorTime               = find(numFactor>= max(numFactor)-1, 1, 'last');             
    fitResult                = fig_sigm((1:numTime)/60, fracActNeuron, factorTime, 0, nan);
    timePoints               = fitResult.x;
    ypred                    = fitResult.ypred;
%     plotAlignPerc(timePoints, ypred, 90); 
    plot(timePoints, ypred); 
    box off
    ylabel('Frac active neuron')
    xlabel('Time from 90% (hour)')
end

%% 3a. radius of communities
function Figure_2_d_3a(clusterList, numTime)    
    factorRadius = zeros(numTime, 1);
    [means, grps] = grpstats(clusterList(:, 3),clusterList(:, 1), {'mean', 'gname'});
    grps         = str2double(grps);
    factorRadius(grps) = means;
    timePoints   = (1:numTime)/60;
    hold on    
    fitResult                = fig_sigm(timePoints, factorRadius, numTime, 0, nan);
    ypred                    = fitResult.ypred;
%     plotAlignPerc(timePoints, ypred, 90); 
    plot(timePoints, ypred); 
    box off
    ylabel('Radius factor')
    xlabel('Time from 90% (hour)')
end

%% 3b. size of communities
function Figure_2_d_3b(clusterList, numTime)    
    factorRadius = zeros(numTime, 1);
    [means, grps] = grpstats(clusterList(:, 2),clusterList(:, 1), {'mean', 'gname'});
    grps         = str2double(grps);
    factorRadius(grps) = means;
    timePoints   = (1:numTime)/60;
    hold on    
    fitResult                = fig_sigm(timePoints, factorRadius, numTime, 0, nan);
    ypred                    = fitResult.ypred;
%     plotAlignPerc(timePoints, ypred, 90);
    plot(timePoints, ypred); 
    ylabel('Size factor')
    xlabel('Time from 90% (hour)')
    box off
end

%% 4. actTime vs location
function Figure_2_d_4(halfActTime, neuronXLoc)
    hold on
    fitActTime = linearFit(neuronXLoc, halfActTime);
    plot(neuronXLoc, fitActTime, 'linewid', 1)
    box off
    ylabel('Frac neuron')
    xlabel('Time (hour)')    
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


function fitResult = fig_sigm(x, y, refDropPoint, lower_, upper_)
    removeIndex      = find(y/max(y) < 0.6);
    factorTime       = refDropPoint;
    removeIndex(removeIndex<factorTime) = [];
    if isempty(removeIndex)
        removeStart  = length(x);
    else
        removeStart  = removeIndex(1);
    end
    init_params      = [0, max(y), x(find(y/max(y)>0.5, 1, 'first')), 1];
    [~, fitResult]   = sigm_fit(x(1:removeStart), y(1:removeStart), [lower_, upper_, nan, nan], init_params, false);
    fitResult.x      = x(1:removeStart)';
end

function yfit   = linearFit(x, y)
    valid_index = ~isnan(x) & ~isnan(y);
    p           = robustfit(x(valid_index),y(valid_index));
    yfit        = polyval(p(end:-1:1), x);
end


function h = plotAlignPerc(x, y, percVal)
    yp     = prctile(y, percVal);
    yInd   = find(y >= yp, 1, 'first');
    if isempty(yInd); keyboard(); end
    h      = plot(x - x(yInd), y, '-');
end

function h = plotAlignByTime(x, y, Time)
end