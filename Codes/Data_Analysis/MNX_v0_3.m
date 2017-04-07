%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5.  MNX neurons -- active neurons
% 
%  EV and activation as a function of time
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function MNX_v0_3(nFile)            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>    
    mnx               = [];
    load([tempDatDir, fileName, '.mat'], 'mnx', 'tracks', 'side', 'leafOrder', 'slicedIndex', 'activeNeuronMat'); 
    load([tempDatDir, 'EV_' fileName, '.mat'], 'EVMat')
    load([tempDatDir, 'LONOLoading_' fileName, '.mat'], 'CorrectedLMat');
    
    if isempty(mnx) || sum(~mnx) == 0
        return;
    end
    
    EVMat             = EVMat(~mnx, :);
    side              = side(~mnx);
    activeNeuronMat   = activeNeuronMat(~mnx, :);
    numTime           = size(EVMat, 2); %#ok<*NODEF>
    numNeuron         = size(EVMat, 1);        
    m                 = ceil(sqrt(numNeuron));
    
    % number of factor a neuron belong to
    % This is always equal to one (one can test it from the following code)
    numFactorPerNeuron= nan(numNeuron, numTime);
    
    for nTime         = 1:numTime
        LMat          = CorrectedLMat{nTime};
        LMat          = LMat(~mnx, :);
        numFactorPerNeuron(:, nTime) = sum(LMat>0, 2);
    end

    figure;
    
    for nNeuron       = 1:numNeuron
        subplot(m, m, nNeuron)
        hold on
        
        if side(nNeuron) == 1
            mColor  = [0    0.4470    0.7410];
        else
            mColor  = [0.6350    0.0780    0.1840];
        end
        
        timeIndex = ~isnan(EVMat(nNeuron, :));
        actCurrNeuron   = activeNeuronMat(nNeuron, :);
        
        plot(find(timeIndex)/60, actCurrNeuron(timeIndex), '.k');
%         plot(find(timeIndex)/60, numFactorPerNeuron(nNeuron, timeIndex), '.m');
        plot(find(timeIndex)/60, EVMat(nNeuron, timeIndex), '.', 'color', mColor);
        
        hold off;
        ylim([0 1])
        xlim([0 numTime/60])    
        xlabel('Time (hour)')
        ylabel('EV')
        box off
    end
    
    setPrint(8*m, 6*m, [plotDir, 'SigFitEVMNXNegActiveNeuron_', fileName], 'pdf');
end