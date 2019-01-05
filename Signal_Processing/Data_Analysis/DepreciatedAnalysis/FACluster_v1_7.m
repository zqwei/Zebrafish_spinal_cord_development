%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.  Covariance analysis: Factor analysis -- number of island using non
% LONO methods with selected neurons
% 
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function FACluster_v1_7(nFile)    

    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>  
    load([tempDatDir, 'FALONO_', fileName, '.mat'], 'EVLONO');
    load([tempDatDir, fileName, '.mat'], 'dff', 'timePoints','activeNeuronMat', 'side');
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
    figure;
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
    
    
    % half time active neuron
    numActNeuron  = mean(activeNeuronMat, 1);
    halfActTimes  = find(numActNeuron>0.49 & numActNeuron<0.51);
    plot(halfActTimes/60, opt1Dim(halfActTimes),'ro')
    
    ylabel('# factors')
    xlabel('Time (hour)')
    set(gca, 'xtick', 1:6)
    xlim([0 numTime/60])
    ylim([0 8])
    box off

    % phase-lock transition time
    activeNeuronMat   = sum(activeNeuronMat,1);
    load([tempDatDir, 'EvoLoading_' fileName, '.mat'], 'networkMat')
    numTime           = length(networkMat); %#ok<USENS>
    networkSummary    = cell(numTime, 1);
    
    for nTime         = 1:numTime
        dat           = [];
        if activeNeuronMat(nTime) > 1
            factorSet     = networkMat{nTime};
            delayMat      = factorSet{end}.delayMat;
            corrMat       = factorSet{end}.corrMat;
            sizeMat       = ones(length(factorSet)-2);
            sideMat       = ones(length(factorSet)-2);
            for nFactor   = 2:length(factorSet)-1
                sizeMat(nFactor-1) = length(factorSet{nFactor}.neuronIndex);
                sideMat(nFactor-1) = mode(side(factorSet{nFactor}.neuronIndex));
            end
            for nFactor   = 1:length(corrMat)
                for mFactor = nFactor+1:length(corrMat)
                    if sizeMat(nFactor)>1 && sizeMat(mFactor)>1 % FA-FA
                        if sideMat(nFactor) ~= sideMat(mFactor) 
                            dat = [dat; delayMat(nFactor, mFactor), corrMat(nFactor, mFactor)];
                        end
                    end
                end
            end
        end
        networkSummary{nTime} = dat;
    end
    
    
    randomMean        = nan(numTime, 1);  
    
    for nTime         = 1:numTime
        dat           = networkSummary{nTime};
        if size(dat, 1)>0
            randomMean(nTime) = mean(isnan(dat(:,1)));
        end
    end        
    
    phaseLockTime = find(randomMean>0, 1, 'last') + 1;
    plot(phaseLockTime/60, opt1Dim(phaseLockTime),'bs')

    setPrint(8, 6, [plotDir, 'numFactorLONOActiveNeuronsTemplate_', fileName], 'pdf');
    
    close all
        
end

function y = doubleSizedGauss(p, x)
    p      =  num2cell(p);
    
    [a, b, c, cr] = deal(p{:});
    d = 0;
    dr = 2;
    
    y = x;
    y(x<=b) = a*exp(-((x(x<=b)-b)/c).^2) + d;
    y(x>b)  = (a+d-dr)*exp(-((x(x>b)-b)/cr).^2) + dr;
end