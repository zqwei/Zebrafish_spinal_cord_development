%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.  Size of cluster at half active time and half EV time
% 
% Half EV - Half active time
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function FAEV_v0_5(nFile)            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>    
    load([tempDatDir, fileName, '.mat'], 'mnx', 'activeNeuronMat'); 
    load([tempDatDir, 'EV_' fileName, '.mat'], 'halfActTime', 'RSquare', 'halfEVTime', 'validFitIndex')
    load([tempDatDir, 'LONOLoading_' fileName, '.mat'], 'CorrectedLMat') 
    
    validNeuronIndex  = (RSquare>0.6 & validFitIndex);
    halfActTimeInt    = ceil(halfActTime*60);
    halfEVTimeInt     = ceil(halfEVTime*60);
    halfActSize       = nan(size(halfActTime));
    halfEVSize        = nan(size(halfEVTime));
    
    
    for nNeuron       = find(validNeuronIndex)'
        halfAct       = halfActTimeInt(nNeuron);
        halfEV        = halfEVTimeInt(nNeuron);
        
        if halfAct>length(CorrectedLMat); halfAct = length(CorrectedLMat); end
        if halfEV>length(CorrectedLMat);  halfEV = length(CorrectedLMat); end
        
        if halfAct<=0; halfAct = 1; end
        if halfEV<=0;  halfEV = 1; end
        
        LMat          = CorrectedLMat{halfAct};
        LMatSize      = sum(LMat>0, 1);
        if sum(LMat(nNeuron, :)>0)>0
            halfActSize(nNeuron) = mean(LMatSize(LMat(nNeuron, :)>0));
        else
            halfActSize(nNeuron) = 1;
        end
        
        
        LMat          = CorrectedLMat{halfEV};
        LMatSize      = sum(LMat>0, 1);
        if sum(LMat(nNeuron, :)>0)>0
            halfEVSize(nNeuron)  = mean(LMatSize(LMat(nNeuron, :)>0));
        else
            halfEVSize(nNeuron)  = 1;
        end
        
    end
    
    
    
    figure;
    hold on
    if ~exist('mnx', 'var')
        plot(halfActTime(RSquare>0.6 & validFitIndex), halfActSize(RSquare>0.6 & validFitIndex), 'ok')
    else
        plot(halfActTime(RSquare>0.6 & validFitIndex & mnx==1), halfActSize(RSquare>0.6 & validFitIndex & mnx==1), 'ok')
        plot(halfActTime(RSquare>0.6 & validFitIndex & mnx==0), halfActSize(RSquare>0.6 & validFitIndex & mnx==0), 'sr')
    end
%     h = refline(1);
%     h.LineStyle = '--';
%     h.Color = [0.5 0.5 0.5];
    hold off
    xlabel('Half active time (hr)')
    ylabel('Size of joint cluster')
    
    setPrint(8, 6, [plotDir, 'HalfActiveTimesFactorSize_', fileName], 'pdf');
    
    
    figure;
    hold on
    if ~exist('mnx', 'var')
        plot(halfEVTime(RSquare>0.6 & validFitIndex), halfEVSize(RSquare>0.6 & validFitIndex), 'ok')
    else
        plot(halfEVTime(RSquare>0.6 & validFitIndex & mnx==1), halfEVSize(RSquare>0.6 & validFitIndex & mnx==1), 'ok')
        plot(halfEVTime(RSquare>0.6 & validFitIndex & mnx==0), halfEVSize(RSquare>0.6 & validFitIndex & mnx==0), 'sr')
    end
%     h = refline(1);
%     h.LineStyle = '--';
%     h.Color = [0.5 0.5 0.5];
    hold off
    xlabel('Half EV time (hr)')
    ylabel('Size of joint cluster')
    
    setPrint(8, 6, [plotDir, 'HalfEVTimesFactorSize_', fileName], 'pdf');
    
    
    close all
end