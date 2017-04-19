% % function Data_Analysis_List_all
numFile = 17;

for nFile = 15 %1:numFile              
% % %     disp(['Generate the orginal data file for dataset #' num2str(nFile)]);   
% % %     disp('==================================');
% % %     disp('Analysis #1 -- clustering analysis');   
% % %     disp('Generate clusters based on correlation matrix'); 
% % %     disp('Plot correlation matrix after being clustered');
% % %     disp('Check the colinearity after being clustered');
% % %     disp('Remove the time points with abrupt change of colinearity');
% % %     Neuron_selection_v0(nFile);
% % %     Neuron_selection_v1(nFile); % Correlation matrix using AP information    
% %     Neuron_selection_v2(nFile); % add New x, y, z coordinates to analysis
% % % 
% % %     disp('==================================');
% % %     disp('Analysis #2 -- factor analysis -- number of factors');   
% % %     disp('Analysis #2 -- factor analysis -- number of factors -- non-LONO'); 
% % %     disp('Analysis #2 -- factor analysis -- number of factors -- LONO'); 
% % %     FACluster_v0_0(nFile) % # factors with selected active neurons
% % %     FACluster_v0_1(nFile) % plot LONO - time-dim (pcolor & contourf)
% % %     FACluster_v0_2(nFile) % plot LONOM optimal number -- LONO or non-LONO
% % %     FACluster_v0_2_1(nFile) % plot LONOM optimal number with Confidence Interval
% % %     FACluster_v0_3(nFile) % compute loading matrix for optimal number FA
% % %     FACluster_v0_4(nFile) % correct FA, removing the projection without low EV
% % %     FACluster_v0_5_1(nFile) % generate FA evolution video -- Yinan version
% %     FACluster_v0_7(nFile) % generate networkMat (delay and correlation mat) for 0_8 and 1_2, 1_3 plots 
% %     FACluster_v0_7_1(nFile) % generate networkMat -- the same as FACluster_v0_7 but for different use of data format
% %     FACluster_v0_9(nFile) % plot FA center and size as a function of time
% % %     FACluster_v1_0(nFile) % generate FA evolution video with convex-hull boundary
% %     FACluster_v1_1(nFile) % plot FA size against time
% %     FACluster_v1_2(nFile) % plot delay time -- FA intra Neuron
% %     FACluster_v1_3(nFile) % plot delay time -- other
% %     FACluster_v1_5(nFile) % plot randomness of contra FA-FA delay time
% %     FACluster_v1_8(nFile) % plot FA-FA delay time without std and mean
% %     FACluster_v1_6(nFile) % plot size and radius of FA as time
% %     
% % 
% % %     disp('==================================');
% % %     disp('Analysis #4 -- factor analysis -- EV for single neurons'); %based on LONO
% % %     FAEV_v0_0(nFile) % generated EV time mat file
% % %     FAEV_v0_1(nFile) % compute half EV time
% % %     FAEV_v0_2(nFile) % compare half EV time vs activation time
% % %     FAEV_v0_3(nFile) % EV time with location
% % 
% % %     disp('==================================');
% % %     disp('Analysis #6 -- MNX');
% %     MNX_v0_0(nFile) % plot num neurons of time as a function of cell type
% %     MNX_v0_1(nFile) % plot half EV time as a function of cell type and location
% %     MNX_v0_2(nFile) % plot half EV time distribution
% %     MNX_v0_3(nFile)
% %     MNX_v0_4(nFile)
% %     MNX_v0_5(nFile)
% %     MNX_v0_6(nFile)
    MNX_v0_7(nFile) % location 
% %     MNX_v0_8(nFile) % contra phase evolutions video
% % %     disp('==================================');
% % %     disp('Analysis #7 -- Phase');
% % %     PhaseCluster_v0_0(nFile)

% % %     close all;
end

