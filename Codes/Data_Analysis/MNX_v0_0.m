%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5.  MNX neurons -- active neurons
% 
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function MNX_v0_0(nFile)            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>    
    load([tempDatDir, fileName, '.mat'], 'mnx', 'activeNeuronMat'); 
    load([tempDatDir, 'LONOLoading_' fileName, '.mat'], 'CorrectedLMat') 
    
    numTime           = size(activeNeuronMat, 2); %#ok<NODEF>
    mnxActiveNeuron   = activeNeuronMat(mnx==0, :);
    sigNeuronsMat     = activeNeuronMat;
    
    for nTime         = 1:numTime
        lambda        = CorrectedLMat{nTime};
        lambda(isnan(lambda)) = 0;
        sigNeuronsMat(:, nTime) = sum(lambda~=0, 2)>0; %#ok<USENS>
    end
    
    mnxSigNeuron      = sigNeuronsMat(mnx==0, :);
    
    figure; 
    hold on
    plot((1:numTime)/60, sum(mnxActiveNeuron, 1), 'ok')
    plot((1:numTime)/60, sum(mnxSigNeuron, 1), 'sr')
    legend({'active', 'factored'})
    legend('location', 'best')
    box off
    xlabel('Time (hour)')
    ylabel('# neurons')
    xlim([0 numTime/60])
%     ylim([0 15])
    setPrint(8, 6, [plotDir, 'MNX-NumNeurons_', fileName], 'pdf');
    
    
    numTime           = size(activeNeuronMat, 2); %#ok<NODEF>
    mnxActiveNeuron   = activeNeuronMat(mnx==1, :);
    sigNeuronsMat     = activeNeuronMat;
    
    for nTime         = 1:numTime
        lambda        = CorrectedLMat{nTime};
        lambda(isnan(lambda)) = 0;
        sigNeuronsMat(:, nTime) = sum(lambda~=0, 2)>0; %#ok<USENS>
    end
    
    mnxSigNeuron      = sigNeuronsMat(mnx==1, :);
    
    figure; 
    hold on
    plot((1:numTime)/60, sum(mnxActiveNeuron, 1), 'ok')
    plot((1:numTime)/60, sum(mnxSigNeuron, 1), 'sr')
    legend({'active', 'factored'})
    legend('location', 'best')
    box off
    xlabel('Time (hour)')
    ylabel('# neurons')
    xlim([0 numTime/60])
%     ylim([0 15])
    setPrint(8, 6, [plotDir, 'MNX+NumNeurons_', fileName], 'pdf');
    
end