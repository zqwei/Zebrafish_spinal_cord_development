%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.  Covariance analysis: Factor analysis -- loadin matrix evolution --
% reordering
% 
% center of the factor
% connectivity strength to each node
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function FACluster_v1_6(nFile) 
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile};     
    load([tempDatDir, fileName, '.mat'], 'activeNeuronMat'); 
    load([tempDatDir, 'EvoLoading_' fileName, '.mat'], 'networkMat', 'neuronXLoc', 'neuronYLoc', 'neuronZLoc'); 
    numTime           = length(networkMat); %#ok<*USENS>
    
    clusterList       = []; % time, size, radius
    
    for nTime         = 1:numTime
        factorSet     = networkMat{nTime}; 
        for nFactor   = 2:length(factorSet)-1
            neuronIndex = factorSet{nFactor}.neuronIndex;
            if length(neuronIndex)>1
                cluster = [ nTime,...
                            length(neuronIndex)];
                xLoc    = neuronXLoc(neuronIndex, nTime);
                yLoc    = neuronYLoc(neuronIndex, nTime);
                zLoc    = neuronZLoc(neuronIndex, nTime);
                radius  = (xLoc - factorSet{nFactor}.x).^2 + (yLoc - factorSet{nFactor}.y).^2 + ...
                          (zLoc - factorSet{nFactor}.z).^2;
                radius  = mean(sqrt(radius));
                cluster = [cluster, radius];
                clusterList = [clusterList; cluster];
            end
        end
    end

    figure;
    plot(clusterList(:, 1)/60, clusterList(:, 3), 'ok', 'markerfacecolor','k')
    box off
    xlabel('Time (hr)')
    ylabel('Radius (um)')
    xlim([0 numTime/60])
    setPrint(8, 6, [plotDir, 'FCRadiusTimeScatter_', fileName], 'pdf')
    
    figure;
    [means, stds, grps] = grpstats(clusterList(:, 3),clusterList(:, 1), {'mean', 'std', 'gname'});
    errorbar(str2double(grps)/60, means, stds, 'ok')
    box off
    xlabel('Time (hr)')
    ylabel('Radius (um)')
    xlim([0 numTime/60])
    setPrint(8, 6, [plotDir, 'FCRadiusTime_', fileName], 'pdf')
    
    figure;
    [means, stds, grps] = grpstats(clusterList(:, 2),clusterList(:, 1), {'mean', 'std', 'gname'});
    errorbar(str2double(grps)/60, means, stds, 'ok')
    box off
    xlabel('Time (hr)')
    ylabel('# Neuron')
    xlim([0 numTime/60])
    setPrint(8, 6, [plotDir, 'FCSizeTime_', fileName], 'pdf')

end