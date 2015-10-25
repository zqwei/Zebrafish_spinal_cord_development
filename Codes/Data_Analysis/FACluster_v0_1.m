%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.  Covariance analysis: Factor analysis -- number of island using non
% LONO methods with all neurons
%
% -- selecting the active neurons...
% 
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function FACluster_v0_1(nFile)            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>    
    load([tempDatDir, fileName, '.mat'], 'dff', 'timePoints');
    numPlot           = length(timePoints);
    numNeuron         = size(dff, 1); %#ok<NODEF>
    
    pValue            = 0.05;
    corrThres         = 0.3;
    sigNeuronsMat     = false(numNeuron, numPlot);
    activeNeuronMat   = false(numNeuron, numPlot);
    
    for nPlot         = 1:numPlot
        slicedDFF     = dff(:,timePoints(nPlot)+1:timePoints(nPlot)+1200);
        slicedDFF     = bsxfun(@minus, slicedDFF, mean(slicedDFF,2));
        slicedDFF     = bsxfun(@rdivide, slicedDFF, std(slicedDFF,[],2))';
        [corrSlicedDFF, sigSlicedDFF] = corr(slicedDFF); % 'type', 'Spearman'        
        corrSlicedDFF(sigSlicedDFF>pValue) = 0;
        sigNeurons    = sum(corrSlicedDFF>corrThres)>0;   
        ksVec         = false(numNeuron, 1); 
        
        for nNeuron   = find(sigNeurons)
            slicedDFFNeuron = slicedDFF(:, nNeuron);
            ksVec(nNeuron) = kstest(slicedDFFNeuron); 
        end     
        
        for nNeuron   = 1:numNeuron
            slicedDFFNeuron = slicedDFF(:, nNeuron);
%             activeNeuronMat(nNeuron, nPlot) = kstest2(...
%                 -slicedDFFNeuron(slicedDFFNeuron<0), slicedDFFNeuron(slicedDFFNeuron>0),...
%                 'alpha', 0.05); %, 'tail', 'smaller'
            activeNeuronMat(nNeuron, nPlot) = kstest(slicedDFFNeuron);
        end
                
        if sum(ksVec) > 1
            slicedDFF     = slicedDFF(:, ksVec);
            [corrSlicedDFF, sigSlicedDFF] = corr(slicedDFF);
            corrSlicedDFF(sigSlicedDFF>pValue) = 0;
            sigNeuronsSub = sum(corrSlicedDFF>corrThres)>0; 
            ksVecIndex    = find(ksVec);
            ksVec(ksVecIndex(~sigNeuronsSub)) = false;
        end
        
        if sum(ksVec) > 1
            sigNeuronsMat(:, nPlot) = ksVec;
        end
    end
    
    save([tempDatDir, fileName, '_numFactorNoLONOActiveNeuronsSimplified.mat'], 'activeNeuronMat', 'sigNeuronsMat', '-append'); 
    
    figure;
    plot(timePoints/4/3600, sum(sigNeuronsMat, 1), 'o')
    xlabel('Time (hour)')
    ylabel('# of active neurons')
    xlim([0 timePoints(end)/4/3600])
    setPrint(8, 6, [plotDir, 'numSigNeurons_', fileName], 'pdf');
    
    figure;
    plot(timePoints/4/3600, sum(activeNeuronMat, 1), 'o')
    xlabel('Time (hour)')
    ylabel('# of active neurons')
    xlim([0 timePoints(end)/4/3600])
    setPrint(8, 6, [plotDir, 'numActNeurons_', fileName], 'pdf');
    
    figure;
    hold on
    plot(timePoints/4/3600, sum(sigNeuronsMat, 1), 's')
    plot(timePoints/4/3600, sum(activeNeuronMat, 1), 'o')
    hold off
    xlabel('Time (hour)')
    ylabel('# of active neurons')
    xlim([0 timePoints(end)/4/3600])
    setPrint(8, 6, [plotDir, 'numActiveNeurons_', fileName], 'pdf');
    
end