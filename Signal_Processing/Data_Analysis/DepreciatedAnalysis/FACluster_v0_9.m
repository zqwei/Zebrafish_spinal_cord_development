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


function FACluster_v0_9(nFile) 
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>    
    load([tempDatDir, fileName, '.mat'], 'activeNeuronMat'); 
    load([tempDatDir, 'EvoLoading_' fileName, '.mat'], 'networkMat', 'neuronXLoc', 'neuronYLoc', 'neuronZLoc'); 
    xtracks           = mean(neuronXLoc);
    ytracks           = mean(neuronYLoc);
    ztracks           = mean(neuronZLoc);
    numTime           = length(xtracks);
    
    clusterList       = []; % location, time, size
    
    for nTime         = 1:numTime
        factorSet     = networkMat{nTime}; 
        x             = xtracks(nTime);
        y             = ytracks(nTime);
        z             = ztracks(nTime);
        for nFactor   = 2:length(factorSet)-1
            neuronIndex = factorSet{nFactor}.neuronIndex;
            if length(neuronIndex)>1
                cluster = [ factorSet{nFactor}.x - x, ...
                            factorSet{nFactor}.y - y, ... 
                            factorSet{nFactor}.z - z, ...
                            nTime,...
                            length(neuronIndex)];
                clusterList = [clusterList; cluster];
            end
        end
    end
       
    figure;
    scatter(clusterList(:, 1), clusterList(:, 2), clusterList(:, 5)*4, clusterList(:, 4), 'filled');
    alpha(0.8)
    colormap(parula)
    h = colorbar;
    ylabel(h, 'Time (min)')
    xlabel('x-axis')
    ylabel('y-axis')
    setPrint(8, 6, [plotDir, 'ClusterCenterXY_', fileName], 'pdf')    
    close all
    
    figure;
    scatter(clusterList(:, 1), clusterList(:, 3), clusterList(:, 5)*4, clusterList(:, 4), 'filled');
    alpha(0.8)
    colormap(parula)
    h = colorbar;
    ylabel(h, 'Time (min)')
    xlabel('x-axis')
    ylabel('z-axis')
    setPrint(8, 6, [plotDir, 'ClusterCenterXZ_', fileName], 'pdf')
    close all

end