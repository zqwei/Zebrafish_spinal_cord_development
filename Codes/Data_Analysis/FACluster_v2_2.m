%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.  Covariance analysis: Factor analysis -- Loading matrix using LONO 
% number of factors -- Varimax -- correction of loading matrix -- comparing
% number of sig. Neurons vs neurons in factors (they should be the same)
% 
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function FACluster_v2_2(nFile)            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>    
    load([tempDatDir, fileName, '.mat'], 'timePoints'); 
    load([tempDatDir, fileName, '_numFactorNoLONOActiveNeuronsSimplified.mat'], 'sigNeuronsMat', 'LONOM')
    load([tempDatDir, fileName, '_LONOLoading.mat'], 'CorrectedLMat');
    
    numTime           = length(LONOM);
    timePoints        = 1:numTime;
    numPlot           = length(timePoints); %#ok<USENS>
    numDiffNeurons    = zeros(numPlot, 1);
    
    for nPlot         = 1:numPlot    
        if sum(sigNeuronsMat(:, nPlot))>1
            LMatnTime     = CorrectedLMat{nPlot};
            numFANeurons  = sum(sum(LMatnTime>0, 2)>0);
            numDiffNeurons(nPlot) = sum(sigNeuronsMat(:, nPlot)) - numFANeurons;
        end
    end

    figure;
    hold on
    plot(timePoints/60, LONOM,'o')
    plot(find(numDiffNeurons>1)/60, LONOM(numDiffNeurons>1),'x')
%     text(find(numDiffNeurons>0)/60, LONOM(numDiffNeurons>0)+0.3, num2str(numDiffNeurons(numDiffNeurons>0)))
    ylabel('# factors')
    xlabel('Time (hour)')
    set(gca, 'xtick', 1:6)
    xlim([0 numTime/60])
    title('Leave one neuron out')
    setPrint(8, 6, [plotDir, 'numFactorLONOMissingActiveNeurons_', fileName], 'pdf');
    hold off
    
    save([tempDatDir, fileName, '_LONOLoading.mat'], 'numDiffNeurons', '-append');
    
end