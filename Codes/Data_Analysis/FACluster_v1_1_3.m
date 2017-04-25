%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.  Covariance analysis: Factor analysis -- loadin matrix evolution --
% reordering
% 
% size of factor against time
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function FACluster_v1_1_3(nFile)            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>    
    load([tempDatDir, fileName, '.mat'], 'activeNeuronMat', 'mnx'); 
    load([tempDatDir, 'EvoLoading_' fileName, '_v2.mat'], 'networkMat')
    numActNeuron      = sum(activeNeuronMat, 1);
    numTime           = length(numActNeuron);
    numFactor         = zeros(numTime, 1);
    factorNeuronMat   = false(size(activeNeuronMat));
    

    
    for nTime         = 1:numTime
        factorSet     = networkMat{nTime, 1};
        for nFactor   = 2:length(factorSet)
            if length(factorSet(nFactor).neuronIndex) > 1
                factorNeuronMat(factorSet(nFactor).neuronIndex, nTime) = true;
                numFactor(nTime)       = numFactor(nTime)+1;
            end
        end
    end
    
    numFactorNeuron   = sum(factorNeuronMat, 1);
            
    figure;
    
    fracNeuron        = numFactorNeuron./numActNeuron;
    
    
    hold on
    plot((1:numTime)/60, fracNeuron, 'ok')
    % remove drop points
    removeIndex              = find(fracNeuron < 0.6);
    factorTime               = find(numFactor>= max(numFactor)-1, 1, 'last');
    removeIndex(removeIndex<factorTime) = [];
    if isempty(removeIndex)
        removeStart          = numTime;
    else
        removeStart          = removeIndex(1);
    end
    % for nfile = 12: removeStart  = 180; fix_fit_max = 0.95;
    % for nfile = 13; fix_fit_max = 0.95;
    % for nfile = 17: removeStart  = 180; fix_fit_max = 0.95;
    % for nfile = 18; removeStart  = 120;
    % for nfile = 19; removeStart  = 150; fix_fit_max = 0.95;
    removeStart  = 180;
    init_params              = [0, max(fracNeuron), find(fracNeuron>0.5, 1, 'first')/60, 1];
    [fitParams, fitResult]   = sigm_fit((1:removeStart)/60, fracNeuron(1:removeStart), [0, 1.0, nan, nan], init_params, false);
    ypred                    = fitResult.ypred;
    ypredlowerCI             = fitResult.ypredlowerCI;
    ypredupperCI             = fitResult.ypredupperCI;
    mColor                   = 'b';
    plot((1:removeStart)/60, ypred, '-', 'linewid', 2.0, 'Color', mColor);
    plot((1:removeStart)/60, ypredlowerCI, '-', 'linewid', 0.5, 'Color', mColor); 
    plot((1:removeStart)/60, ypredupperCI, '-', 'linewid', 0.5, 'Color', mColor); 
    ylim([0 1])
    gridxy([fitParams(3)-1/fitParams(4)*log10(9) fitParams(3)+1/fitParams(4)*log10(9)], [], 'color', 'r', 'linestyle', '--')
    hold off
    xlabel('Time (hour)')
    ylabel('frac. neuron factored')
    xlim([0 numTime/60])
    
    box off
    disp([num2str(nFile) ':' fileName ': ' num2str(1/fitParams(4)*(log10(9)+log10(9)))])
    
    setPrint(8, 6, [plotDir, 'FASizeSigmoidFit_', fileName], 'pdf');
        
end