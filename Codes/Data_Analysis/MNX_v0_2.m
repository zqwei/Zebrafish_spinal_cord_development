%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5.  MNX neurons -- active neurons
% 
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function MNX_v0_2(nFile)            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>    
    load([tempDatDir, fileName, '.mat'], 'mnx', 'side','tracks'); 
    load([tempDatDir, 'EV_' fileName, '.mat'], 'halfEVTime', 'RSquare', 'halfActTime', 'validFitIndex')
    
    RSquareThres      = 0.6;
    
    maxTime           = ceil(max(halfEVTime(RSquare > RSquareThres & validFitIndex)));
    
    
    figure;
    hold on
    [f, xi] = ksdensity(halfEVTime(mnx==0 & RSquare > RSquareThres & validFitIndex), 0:0.1:maxTime);
    plot(xi, f)
    [f, xi] = ksdensity(halfEVTime(mnx==1 & RSquare > RSquareThres & validFitIndex), 0:0.1:maxTime);
    plot(xi, f)
    legend({'MNX-', 'MNX+'})
    legend('location', 'best')
    legend('boxoff')
    xlabel('Half EV time (hour)')
    ylabel('Prob. Density');
    setPrint(8, 6, [plotDir, 'MNXHalfEVDistr_', fileName], 'pdf');
    close all

    
    figure;
    maxTime           = ceil(max(halfActTime(RSquare > RSquareThres & validFitIndex)));
    hold on
    [f, xi] = ksdensity(halfActTime(mnx==0 & RSquare > RSquareThres & validFitIndex), 0:0.1:maxTime);
    plot(xi, f)
    [f, xi] = ksdensity(halfActTime(mnx==1 & RSquare > RSquareThres & validFitIndex), 0:0.1:maxTime);
    plot(xi, f)
    legend({'MNX-', 'MNX+'})
    legend('location', 'best')
    legend('boxoff')
    xlabel('Half Activation time (hour)')
    ylabel('Prob. Density');
    setPrint(8, 6, [plotDir, 'MNXHalfActDistr_', fileName], 'pdf');
    close all
    
    figure;
    halfDifferenceTime = halfEVTime - halfActTime;
    maxTime           = ceil(max(halfDifferenceTime(RSquare > RSquareThres & validFitIndex)));
    minTime           = floor(min(halfDifferenceTime(RSquare > RSquareThres & validFitIndex)));
    hold on
    [f, xi] = ksdensity(halfDifferenceTime(mnx==0 & RSquare > RSquareThres & validFitIndex), minTime:0.1:maxTime);
    plot(xi, f)
    [f, xi] = ksdensity(halfDifferenceTime(mnx==1 & RSquare > RSquareThres & validFitIndex), minTime:0.1:maxTime);
    plot(xi, f)
    legend({'MNX-', 'MNX+'})
    legend('location', 'best')
    legend('boxoff')
    xlabel('t_{EV} - t_{act} (hour)')
    ylabel('Prob. Density');
    setPrint(8, 6, [plotDir, 'MNXHalfDiffTimeDistr_', fileName], 'pdf');
    close all

    
end