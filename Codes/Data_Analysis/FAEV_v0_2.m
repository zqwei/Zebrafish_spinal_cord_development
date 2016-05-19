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
    load([tempDatDir, 'EV_' fileName, '.mat'], 'halfActTime', 'RSquare', 'halfEVTime', 'validFitIndex')
    
    figure;
    hold on
    plot(halfActTime(RSquare>0.6 & validFitIndex), halfEVTime(RSquare>0.6 & validFitIndex), 'ok')
    h = refline(1);
    h.LineStyle = '--';
    h.Color = [0.5 0.5 0.5];
    hold off
    xlabel('Half active time (hr)')
    ylabel('Half EV time (hr)')
    
    setPrint(8, 6, [plotDir, 'HalfEVActiveNeuronTimes_', fileName]);
    
end