% % function Data_Analysis_List_Ablation
addpath('../Data_Analysis/')
fishList = [19,22,11,25,26,20,24,23];


% LMatThres = 0.3;
% PsiMatThres = 0.8;
% LifetimeThres = 0;
% for i = 1:numel(fishList)
%     nFile = 24+fishList(i)*2-1;
%     disp(['Generate the orginal data file for dataset #' num2str(nFile)]);
% %     Neuron_selection_v0_short_win(nFile);
% % 
% %     Neuron_selection_v1(nFile); % Correlation matrix using AP information
% %     Neuron_selection_v2(nFile); 
% % 
% %     FACluster_v0_0(nFile) % # factors with selected active neurons
% %     FACluster_v0_2_short_win(nFile) % plot LONOM optimal number -- LONO or non-LONO
% % 
% %     FACluster_v0_3_short_win(nFile) % compute loading matrix for optimal number FA
% %     FACluster_v0_4(nFile) % correct FA, removing the projection without low EV
% %     
% %     FACluster_v0_3_whole_win(nFile)
% %     FACluster_v0_5_1_short_win(nFile) % generate FA evolution video
%     FACluster_v0_3_stable_win(nFile, LMatThres, PsiMatThres, LifetimeThres);
% 
% 
% %     close all;
% end

% Ablation_v_4_2(fishList); %plot change/ratio of synchronization level

Ablation_v_5_1(fishList); % plot experimental summary in matrix