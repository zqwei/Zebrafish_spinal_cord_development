% Fit for Figure 2c for each fish
% 1. number of communities
% 2. Fraction of non-factored neurons
% 2. percentage of total active neurons
% 3. radius of communities
% 3. size of communities


function Figure_2_d(nFile)


end




function Figure_2_d_1(nFile)    

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