% % function Data_Analysis_List_all
% numFile = 19;

fileToAnalysis = [3 4 5 7 10 11 12 13 15 16 17 18 19];

for nFile = fileToAnalysis           
% % %     Neuron_selection_v0(nFile);
% % %     Neuron_selection_v1(nFile); % Correlation matrix using AP information    
% % %     Neuron_selection_v2(nFile); % add New x, y, z coordinates to analysis

% % %     FACluster_v0_0(nFile) % # factors with selected active neurons
% % %     FACluster_v0_1(nFile) % plot LONO - time-dim (pcolor & contourf)
% % %     FACluster_v0_2(nFile) % plot LONOM optimal number -- LONO or non-LONO
% % %     FACluster_v0_2_1(nFile) % plot LONOM optimal number with Confidence Interval
% % %     FACluster_v0_3(nFile) % compute loading matrix for optimal number FA
% % %     FACluster_v0_4(nFile) % correct FA, removing the projection without low EV
% % %     FACluster_v0_2_2(nFile) % plot LONOM optimal number with dropping the overlapped factors
% % %     FACluster_v0_5_1(nFile) % generate FA evolution video -- Yinan version
% % %     FACluster_v0_7(nFile) % generate networkMat (delay and correlation mat) for 0_8 and 1_2, 1_3 plots 
% % %     FACluster_v0_7_1(nFile) % generate networkMat -- the same as FACluster_v0_7 but for different use of data format
% % %     FACluster_v0_9_1(nFile) % plot FA center and size as a function of time : new coordinates
% % %     FACluster_v1_1(nFile) % plot FA size against time
% % %     FACluster_v1_1_1(nFile) % plot max FA size for each side against time
% % %     FACluster_v1_1_2(nFile) % plot FA size with mnx factored time
    FACluster_v1_1_3(nFile) % plot FA size with mnx factored time
% % %     FACluster_v1_2(nFile) % plot delay time -- FA intra Neuron
% % %     FACluster_v1_3(nFile) % plot delay time -- other
% % %     FACluster_v1_5(nFile) % plot randomness of contra FA-FA delay time
% % %     FACluster_v1_8(nFile) % plot FA-FA delay time without std and mean
% % %     FACluster_v1_6(nFile) % plot size and radius of FA as time
    
% % %     FAEV_v0_0(nFile) % generated EV time mat file
% % %     FAEV_v0_1(nFile) % compute half EV time
% % %     FAEV_v0_2(nFile) % compare half EV time vs activation time
% % %     FAEV_v0_3(nFile) % EV time with location

% % %     MNX_v0_0(nFile) % plot num neurons of time as a function of cell type
% % %     MNX_v0_1(nFile) % plot half EV time as a function of cell type and location
% % %     MNX_v0_2(nFile) % plot half EV time distribution
% % %     MNX_v0_3(nFile)
% % %     MNX_v0_4(nFile)
% % %     MNX_v0_6(nFile)
% % %     MNX_v0_7(nFile) % location 
% % %     MNX_v0_8(nFile) % contra phase evolutions video
% % %     MNX_v0_9(nFile) % generate FA evolution video -- Yinan version with mnx- tag

%     PhaseCluster_v0_0(nFile)

    close all;
end

