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


function FACluster_v0_9_1(nFile) 
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>    
    load([tempDatDir, fileName, '.mat'], 'activeNeuronMat'); 
%     load([tempDatDir, 'EvoLoading_' fileName, '.mat'], 'networkMat', 'neuronXLoc', 'neuronYLoc', 'neuronZLoc'); 
    if ~exist([tempDatDir, 'EvoLoading_' fileName, '_v2.mat'], 'file'); return; end
    load([tempDatDir, 'EvoLoading_' fileName, '_v2.mat'], 'networkMat', 'neuronXLoc', 'neuronYLoc', 'neuronZLoc'); 
    xtracks           = mean(neuronXLoc);
    ytracks           = mean(neuronYLoc);
    ztracks           = mean(neuronZLoc);
    z_max             = max(neuronZLoc(:));
    numTime           = length(xtracks);
    
    clusterList       = []; % location, time, size
    
    for nTime         = 1:numTime
        factorSet     = networkMat{nTime, 1};  %#ok<*USENS>
        x             = xtracks(nTime);
        y             = ytracks(nTime);
        z             = ztracks(nTime);
        for nFactor   = 2:length(factorSet)
            neuronIndex = factorSet(nFactor).neuronIndex;
            if length(neuronIndex)>1
                cluster = [ factorSet(nFactor).x, ...
                            factorSet(nFactor).y, ... 
                            factorSet(nFactor).z/z_max*1.8, ...
                            nTime,...
                            length(neuronIndex)];
                clusterList = [clusterList; cluster];
            end
        end
    end
       
    figure;
    subplot(3, 1, 1)
    scatter(clusterList(:, 1), clusterList(:, 2), clusterList(:, 5)*4, clusterList(:, 4), 'filled');
    alpha(0.8)
    colormap(parula)
    xlim([0 10])
    ylim([-2 2])
    gridxy([0:10], [0], 'color', 'k', 'linestyle', '--')
    h = colorbar;
    ylabel(h, 'Time (min)')
    xlabel('x-axis')
    ylabel('y-axis')
    
    subplot(3, 1, 2)
    side_ind = clusterList(:, 2)>0;
    scatter(clusterList(side_ind, 1), clusterList(side_ind, 3), clusterList(side_ind, 5)*4, clusterList(side_ind, 4), 'filled');
    alpha(0.8)
    colormap(parula)
    xlim([0 10])
    ylim([0 1])
    gridxy([0:10], [], 'color', 'k', 'linestyle', '--')
    h = colorbar;
    ylabel(h, 'Time (min)')
    xlabel('x-axis')
    ylabel('left z-axis')

    subplot(3, 1, 3)
    side_ind = clusterList(:, 2)<0;
    scatter(clusterList(side_ind, 1), clusterList(side_ind, 3), clusterList(side_ind, 5)*4, clusterList(side_ind, 4), 'filled');
    alpha(0.8)
    colormap(parula)
    xlim([0 10])
    ylim([0 1])
    gridxy([0:10], [], 'color', 'k', 'linestyle', '--')
    h = colorbar;
    ylabel(h, 'Time (min)')
    xlabel('x-axis')
    ylabel('right z-axis')
    
    
    setPrint(18, 18, [plotDir, 'ClusterCenter_', fileName], 'pdf')    
    
% %     figure;
% %     scatter(clusterList(:, 1), clusterList(:, 3), clusterList(:, 5)*4, clusterList(:, 4), 'filled');
% %     alpha(0.8)
% %     colormap(parula)
% %     h = colorbar;
% %     ylabel(h, 'Time (min)')
% %     xlabel('x-axis')
% %     ylabel('z-axis')
% %     setPrint(8, 6, [plotDir, 'ClusterCenterXZ_', fileName], 'pdf')
% %     close all

end