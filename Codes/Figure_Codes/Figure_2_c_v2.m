% Fit for Figure 2c for each fish
% 1. number of communities
% 2a. Fraction of non-factored neurons
% 2b. percentage of total active neurons
% 3a. radius of communities
% 3b. size of communities
% 4.  (alternative) time vs. ClusterCenterMaxLoc
% 5. mnx+/mnx- vs. actTime/patternTime

function stats = Figure_2_c_v2(nFile, tag)
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile};   
    load([tempDatDir, 'FALONO_', fileName, '.mat'], 'dumpDuplicatedFactorLONOM');
    load([tempDatDir, fileName, '.mat'],'activeNeuronMat', 'side', 'mnx');
    if ~exist([tempDatDir, 'EvoLoading_' fileName, '_v2.mat'], 'file'); return; end
    load([tempDatDir, 'EvoLoading_' fileName, '_v2.mat'], 'networkMat', 'neuronXLoc', 'neuronYLoc', 'neuronZLoc')
    load([tempDatDir, 'Leader_', fileName, '.mat'], 'activeTime', 'patternTime');

    
    numTime           = length(networkMat); %#ok<*USENS>
    clusterList       = []; % time, size, radius, mean xloc  
    maxClusterXLoc    = nan(numTime, 1);
    for nTime         = 1:numTime
        factorSet     = networkMat{nTime, 1}; 
        maxFacSize    = -1;
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
%                 radius  = (xLoc - factorSet(nFactor).x).^2;
%                 radius  = sqrt(mean(radius));
                if length(neuronIndex)>maxFacSize
                    maxClusterXLoc(nTime) = mean(xLoc);
                    maxFacSize            = length(neuronIndex);
                end
                radius = max(xLoc) - min(xLoc);
                cluster = [cluster, radius];
                clusterList = [clusterList; cluster];
            end
        end
    end
    
    totPlots = 6;
    stats = cell(totPlots, 1);
    
    figure('Position', [0, 0, 1200, 200]);
    
    subplot(1, totPlots, 1)
    stats{1} = Figure_2_c_1(dumpDuplicatedFactorLONOM);    
    
    subplot(1, totPlots, 2)
    stats{2} = Figure_2_c_2(activeNeuronMat, networkMat);
    
    subplot(1, totPlots, 3)
    stats{3} = Figure_2_c_3a(clusterList, numTime);
    
    subplot(1, totPlots, 4)
    stats{4} = Figure_2_c_3b(clusterList, numTime);
    
    subplot(1, totPlots, 5)
    stats{5} = Figure_2_c_4(maxClusterXLoc);
    
    if strcmp(tag, 'actTime')
        subplot(1, totPlots, 6)
        stats{6} = Figure_2_c_5(activeTime, mnx);
    elseif strcmp(tag, 'patternTime')
        subplot(1, totPlots, 6)
        stats{6} = Figure_2_c_5(patternTime, mnx);
    end
    
    
    
    
    setPrint(8*totPlots, 6, [plotDir 'Figure_2b_v1_' fileName '_' tag '_hp'], 'pdf')
    
%     close(gcf)
end

%% 1. number of communities
function pred = Figure_2_c_1(LONOM)    
    numTime           = numel(LONOM);
    timePoints        = (1:numTime)';
    hold on    
    plot(timePoints/60, LONOM,'o', 'color', [0.7 0.7 0.7])    
        
%     % option1: double-gaussian on curve
%     % trancate super-long movies, pad 2 to the end of super-short movies
%     fitLONOM = LONOM;
%     if numTime > 200
%         fitLONOM = fitLONOM(1:200);
%     end
%     fitLONOM(end+1:end+200) = 2;
%     fitTimePoints = (1:numel(fitLONOM))';
%     % fit of numFactor curve
%     fitResult    = fit(fitTimePoints/60, fitLONOM-2, 'gauss1');
%     b            = fitResult.b1;
%     a            = fitResult.a1;
%     c            = fitResult.c1;
%     cr           = c;
%     fitResult    = lsqcurvefit(@(p, x) doubleSizedGauss(p, x), [a, b, c, cr], fitTimePoints/60, fitLONOM);    
%     opt1Dim      = doubleSizedGauss(fitResult,fitTimePoints/60);  
    

%   % option 2: pad 0 at beginning and do sgolay smoothing
%     tmpLONOM = [zeros(100, 1); LONOM; ones(100, 1)*2];
%     opt1Dim  = smooth(tmpLONOM, 201, 'sgolay', 5);
%     opt1Dim  = opt1Dim(101:end);

    % option 3: hpfilt
    tmpLONOM = [zeros(100, 1); LONOM; ones(100, 1)*2];
    lambda  = l1tf_lambdamax(tmpLONOM);
    opt1Dim = hpfilter(tmpLONOM, lambda*2);
    opt1Dim  = opt1Dim(101:end);
%     lambda  = l1tf_lambdamax(LONOM);
%     opt1Dim  = hpfilter(LONOM, lambda*2);


    opt1Dim      = opt1Dim(timePoints);
    plot(timePoints/60, opt1Dim,'k-', 'linewid', 2)
    xlim([0 max(timePoints)/60])  
    ylim([-1, 6]);
    ylabel('Num factor')
    xlabel('Time (hour)')
    box off
    pred.t = timePoints/60;
    pred.y = opt1Dim;
    set(gca, 'TickDir', 'out');
end

%% 2. 
% a. Fraction of non-factored neurons 
% b. percentage of total active neurons
function pred = Figure_2_c_2(activeNeuronMat, networkMat)    
    numActNeuron      = sum(activeNeuronMat, 1);
    fracActNeuron     = numActNeuron/sum(sum(activeNeuronMat, 2)>0);    %exclude never-active neurons
%     fracActNeuron     = sum(activeNeuronMat(sum(activeNeuronMat, 2)>20, :), 1)/sum(sum(activeNeuronMat, 2)>40);    %exclude short-lived neurons
    numTime           = length(numActNeuron);
    numFactor         = zeros(numTime, 1);
    factorNeuronMat   = false(size(activeNeuronMat));
    % compute number of factored neuron
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
    plot((1:numTime)/60, fracActNeuron, 'or')
    plot((1:numTime)/60, 1 - fracNeuron, 'ok')
    
    % fit for fracActNeuron    
    factorTime               = find(numFactor>= max(numFactor)-1, 1, 'last');             
    fitResult                = fig_sigm((1:numTime)/60, fracActNeuron, factorTime, 0, nan);
    timePoints               = fitResult.x;
    ypred                    = fitResult.ypred;
    ypredlowerCI             = fitResult.ypredlowerCI;
    ypredupperCI             = fitResult.ypredupperCI;
    plot(timePoints, ypred, '-', 'linewid', 2.0, 'Color', 'r');
    plot(timePoints, ypredlowerCI, '-', 'linewid', 0.5, 'Color', 'r');
    plot(timePoints, ypredupperCI, '-', 'linewid', 0.5, 'Color', 'r');
    ylim([0 1])
    xlim([0 max(timePoints)])    
    pred.t1 = timePoints;
    pred.y1 = ypred;

    
    % fit for fracNeuron           
    fitResult                = fig_sigm((1:numTime)/60, fracNeuron, factorTime, 0, 1);
    timePoints               = fitResult.x;
    ypred                    = fitResult.ypred;
    ypredlowerCI             = fitResult.ypredlowerCI;
    ypredupperCI             = fitResult.ypredupperCI;
    plot(timePoints, 1- ypred, '-', 'linewid', 2.0, 'Color', 'k');
    plot(timePoints, 1- ypredlowerCI, '-', 'linewid', 0.5, 'Color', 'k');
    plot(timePoints, 1- ypredupperCI, '-', 'linewid', 0.5, 'Color', 'k');
    ylim([0 1])
    xlim([0 max(timePoints)])  
    ylabel('Frac neuron')
    xlabel('Time (hour)')    
    box off
    
    pred.t2 = timePoints;
    pred.y2 = 1- ypred;
    pred.fracAct50 = fitResult.param(1);
    set(gca, 'TickDir', 'out');
end

%% 3a. radius of communities
function pred = Figure_2_c_3a(clusterList, numTime)    
    factorRadius = zeros(numTime, 1);
    [means, stds, grps] = grpstats(clusterList(:, 3),clusterList(:, 1), {'mean', 'std', 'gname'});
    grps         = str2double(grps);
    factorRadius(grps) = means;
    factorStd    = zeros(numTime, 1);
    factorStd(grps) = stds; 
    timePoints   = (1:numTime)/60;
    hold on
    plot(timePoints, factorRadius, 'ok')
%     errorbar(timePoints, factorRadius, factorStd, '.k')
    
    fitResult                = fig_sigm(timePoints, factorRadius, numTime, 0, nan);
    ypred                    = fitResult.ypred;
    ypredlowerCI             = fitResult.ypredlowerCI;
    ypredupperCI             = fitResult.ypredupperCI;
    plot(timePoints, ypred, '-', 'linewid', 2.0, 'Color', 'k');
    plot(timePoints, ypredlowerCI, '-', 'linewid', 0.5, 'Color', 'k');
    plot(timePoints, ypredupperCI, '-', 'linewid', 0.5, 'Color', 'k');
    xlim([0 max(timePoints)])  
    ylabel('Radius factor')
    xlabel('Time (hour)')
    box off
    pred.t = timePoints;
    pred.y = ypred;
    set(gca, 'TickDir', 'out');
end

%% 3b. size of communities
function pred = Figure_2_c_3b(clusterList, numTime)    
    factorRadius = zeros(numTime, 1);
    [means, grps] = grpstats(clusterList(:, 2),clusterList(:, 1), {'mean', 'gname'});
    grps         = str2double(grps);
    factorRadius(grps) = means;
    timePoints   = (1:numTime)/60;
    
    fitResult                = fig_sigm(timePoints, factorRadius, numTime, 0, nan);
    factorRadius = factorRadius/fitResult.param(1); % normalize by sigmoid amplitude
    
    hold on
%     plot(timePoints, factorRadius, 'ok')
    plot(clusterList(:, 1)/60, clusterList(:, 2), 'ok');
    fitResult                = fig_sigm(timePoints, factorRadius, numTime, 0, nan);
    ypred                    = fitResult.ypred;
%     ypredlowerCI             = fitResult.ypredlowerCI;
%     ypredupperCI             = fitResult.ypredupperCI;
%     plot(timePoints, ypred, '-', 'linewid', 2.0, 'Color', 'k');
%     plot(timePoints, ypredlowerCI, '-', 'linewid', 0.5, 'Color', 'k');
%     plot(timePoints, ypredupperCI, '-', 'linewid', 0.5, 'Color', 'k');
    xlim([0 max(timePoints)])  
    ylabel('Factor size')
    xlabel('Time (hour)')
    box off
    pred.t = timePoints;
    pred.y = ypred;
    set(gca, 'TickDir', 'out');
end


% %% 4. actTime/patternTime vs location, boxplot
% function pred = Figure_2_c_4(time, neuronXLoc)
%     neuronXLoc = neuronXLoc(time>0.1 & ~isnan(time));
%     time       = time(time>0.1 & ~isnan( time));
%     segXLoc    = round(neuronXLoc);
%     segNum     = unique(segXLoc);
%     boxplot(time, segXLoc, 'Positions', segNum);
%     hold on
%     
% %     % use group stats
% %     [segTime, segCount]    = grpstats(time, segXLoc, {'mean', 'numel'});
% %     segNum(segCount<=1) = [];
% %     segTime(segCount<=1) = [];
% 
%     
%     % use raw stats
%     segNum  = neuronXLoc;
%     segTime = time;
%     
%     [rho, pval] = corr(segNum, segTime, 'rows', 'pairwise', 'type', 'spearman');
%     fitActTime = linearFit(segNum, segTime);
%     plot(segNum, fitActTime, '-k', 'linewid', 1)
%     box off
%     xlabel('x location (segments)')
%     ylabel('Pattern time (hour)')
%     pred.t = [min(segNum), max(segNum)];
%     pred.y = [min(fitActTime), max(fitActTime)];
%     pred.pval = pval;
%     set(gca, 'TickDir', 'out');
%     ylim([0, ceil(nanmax(time))])
%     title(['rho=' num2str(rho) ', pval=' num2str(pval)]);
% end

%% alternative c_4. time vs. ClusterCenterMaxLoc
function pred = Figure_2_c_4(maxClusterXLoc)
    timePoints = (1:numel(maxClusterXLoc))'/60;
    timePoints(isnan(maxClusterXLoc)) = [];
    maxClusterXLoc(isnan(maxClusterXLoc)) = [];
    
    [rho, pval] = corr(timePoints, maxClusterXLoc, 'type', 'pearson');
    fitClusterXLoc = linearFit(timePoints, maxClusterXLoc);
    hold on
    plot(timePoints, maxClusterXLoc, 'ok')
    plot(timePoints, fitClusterXLoc, '-k', 'linewid', 1)
    box off
    xlabel('time (hour)')
    ylabel('Main factor AP (seg)')
    pred.t = [min(timePoints), max(timePoints)];
    pred.y = [min(fitClusterXLoc), max(fitClusterXLoc)];
    pred.pval = pval;
    set(gca, 'TickDir', 'out');
    xlim([0 max(timePoints)])  
    title(['rho=' num2str(rho) ', pval=' num2str(pval)]);
end

%% 5. mnx+/mnx- vs. actTime/patternTime
function pred = Figure_2_c_5(time, mnx)
    mnx = mnx(~isnan(time));
    time = time(~isnan(time));
    if sum(mnx==0) == 0
        pred.t = [];
        pred.y = [];
        pred.pval = [];
        return;
    end
%     distributionPlot([y1, y2], 'color', [0.5, 0.5, 1], 'showMM', 0, 'globalNorm', 3);
    y1 = time(mnx>0 & ~isnan(time));
    y2 = time(mnx==0 & ~isnan(time));
    
    pval = ranksum(y1, y2);
    
    hold on
    scatter(ones(numel(y1), 1)+randn(numel(y1), 1)/10, y1,10, 'filled', 'MarkerFaceColor', [0.8, 0.8, 0.8]);
    scatter(ones(numel(y2), 1)+randn(numel(y2), 1)/10+1, y2, 10, 'filled', 'MarkerFaceColor', [0.8, 0.8, 0.8]);
    h = boxplot(time, 1-mnx, 'labels', {'mnx+', 'mnx-'});
    set(h, 'linew', 1);
    hold off

    ylabel('Activation time (hour)')
    box off
    set(gca, 'TickDir', 'out');
    title(['pval=' num2str(pval)]);
    pred.t = [1, 2];
    pred.y = [median(y1), median(y2)];
    pred.pval = pval;
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
%     y                = medfilt1(y);
    removeIndex      = find(y/max(y) < 0.6);
    factorTime       = refDropPoint;
    removeIndex(removeIndex<factorTime) = [];
    if isempty(removeIndex)
        removeStart  = length(x);
    else
        removeStart  = removeIndex(1);
    end
    init_params      = [0, max(y), x(find(y/max(y)<0.5, 1, 'last')), 1];
    [~, fitResult]   = sigm_fit(x(1:removeStart), y(1:removeStart), [lower_, upper_, nan, nan], init_params, false);
    fitResult.x      = x(1:removeStart)';
end

function yfit   = linearFit(x, y)
    valid_index = ~isnan(x) & ~isnan(y);
    p           = robustfit(x(valid_index),y(valid_index));
    yfit        = polyval(p(end:-1:1), x);
end