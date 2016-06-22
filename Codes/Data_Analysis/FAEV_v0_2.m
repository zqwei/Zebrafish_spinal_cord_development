%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.  Explained variance: Factor analysis -- corrected loading matrix --
% sigmoid fit
% 
% Half EV - Half active time
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function FAEV_v0_2(nFile)            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>    
    load([tempDatDir, fileName, '.mat'], 'mnx', 'activeNeuronMat'); 
    load([tempDatDir, 'EV_' fileName, '.mat'], 'halfActTime', 'RSquare', 'halfEVTime', 'validFitIndex')
    
    figure;
    hold on
    if ~exist('mnx', 'var')
        plot(halfActTime(RSquare>0.6 & validFitIndex), halfEVTime(RSquare>0.6 & validFitIndex), 'ok')
    else
        plot(halfActTime(RSquare>0.6 & validFitIndex & mnx==1), halfEVTime(RSquare>0.6 & validFitIndex & mnx==1), 'ok')
        plot(halfActTime(RSquare>0.6 & validFitIndex & mnx==0), halfEVTime(RSquare>0.6 & validFitIndex & mnx==0), 'sr')
    end
    h = refline(1);
    h.LineStyle = '--';
    h.Color = [0.5 0.5 0.5];
    hold off
    xlabel('Half active time (hr)')
    ylabel('Half EV time (hr)')
    
    setPrint(8, 6, [plotDir, 'HalfEVActiveNeuronTimes_', fileName]);
    
end