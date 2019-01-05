%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.  Covariance analysis: Factor analysis -- loadin matrix evolution --
% reordering
%
% center of the factor: movement in time along x,y,z
%
%
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
%


function FACluster_v0_9_2(nFile)
addpath('../Func');
setDir;
fileName          = fileNames{nFile};
load([tempDatDir, fileName, '.mat'], 'activeNeuronMat', 'new_x', 'new_y', 'new_z');
load([tempDatDir, 'EvoLoading_' fileName, '_v2.mat'], 'networkMat');
numTime           = size(networkMat, 1);

clusterList       = []; % location, time, size
maxFactorList     = [];
z_max = max(new_z);
x_max = max(new_x);


for nTime         = 1:numTime
    factorSet     = networkMat{nTime, 1};  %#ok<*USENS>
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

c_max = max(clusterList(:,  5));
sizeFactor = 2;
%% visualizing centers with size indicator
figure;
subplot(3, 1, 1);
scatter(clusterList(:, 4)/60, clusterList(:, 1), clusterList(:, 5)*sizeFactor,clusterList(:, 5), 'filled');
ylabel('AP location')
xlabel('Time (h)')
ylim([0, ceil(x_max)]);
caxis([0, c_max]);
hcb = colorbar;
ylabel(hcb, 'Factor size');

subplot(3, 1, 2);
sideTag = clusterList(:, 2)>0;
scatter(clusterList(sideTag, 4)/60, clusterList(sideTag, 3), clusterList(sideTag, 5)*sizeFactor,clusterList(sideTag, 5), 'filled');
ylabel('DV location - left')
xlabel('Time (h)')
caxis([0, c_max]);
ylim([0 1]);
hcb = colorbar;
ylabel(hcb, 'Factor size');
subplot(3, 1, 3);
scatter(clusterList(~sideTag, 4)/60, clusterList(~sideTag, 3), clusterList(~sideTag, 5)*sizeFactor,clusterList(~sideTag, 5), 'filled');
ylabel('DV location - right')
xlabel('Time (h)')
caxis([0, c_max]);
ylim([0, 1]);
hcb = colorbar;
ylabel(hcb, 'Factor size');


setPrint(18, 18, [plotDir, 'ClusterCenterXZ1D_', fileName], 'pdf')

maxFactorList = [];
for nTime = 1:numTime
    cluster = clusterList(clusterList(:, 4)==nTime, :);
    maxFactorIndex = cluster(:, 5)==max(cluster(:, 5));
    xLocMaxFactor = cluster(maxFactorIndex, 1);
    maxFactorList = [maxFactorList; [ones(size(xLocMaxFactor))*nTime/60, xLocMaxFactor]];    
end
figure,
scatter(maxFactorList(:, 1), maxFactorList(:, 2), 'ok');
[rho, pval] = corr(maxFactorList(:, 1), maxFactorList(:, 2), 'type', 'spearman');
title(['rho=', num2str(rho), ',pval=', num2str(pval)]);
setPrint(8, 6, [plotDir, 'ClusterCenterMaxLoc_', fileName], 'pdf')
end