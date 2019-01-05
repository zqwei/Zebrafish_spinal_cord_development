% % function Data_Analysis_List_Ablation
addpath('../Data_Analysis/')
% fishList = [19,22,11,25,26,20,24,23];

% fishList = 27;

fishList = [29,32,22, 30,31,33, 34];
% fishListAll = [22, 29, 30, 31, 32, 33, 34, 11, 18, 20, 23:27, 35:36];
fishListAll = [18,35,29,32,22, 30,31,33, 34, 25, 26, 20, 11, 23, 27, 36, 24];

% %% FA 
% for i = 1:numel(fishList)
%     for nFile = 23+fishList(i)*2 + [0, 1];
%         disp(['Generate the orginal data file for dataset #' num2str(nFile)]);
%         Neuron_selection_v0_short_win(nFile);
% 
%         Neuron_selection_v1(nFile); % Correlation matrix using AP information
%         Neuron_selection_v2(nFile); 
% 
%         FACluster_v0_0(nFile) % # factors with selected active neurons
%         FACluster_v0_2_short_win(nFile) % plot LONOM optimal number -- LONO or non-LONO
% 
%         FACluster_v0_3_short_win(nFile) % compute loading matrix for optimal number FA
%         FACluster_v0_4(nFile) % correct FA, removing the projection without low EV
% 
%     close all;
%     end
% end
% 
% %% Analysis of FA result
% LMatThres = 0.3;
% PsiMatThres = 0.5;
% LifetimeThres = 1; % factor life threshold only applies to after ablation
% for i = 1:numel(fishListAll)
%    for  nFile = 23+fishListAll(i)*2
% %         FACluster_v0_3_whole_win(nFile, LMatThres, PsiMatThres)
% %         FACluster_v0_5_1_short_win(nFile, LMatThres, PsiMatThres) % generate FA evolution video
%         FACluster_v0_3_stable_win(nFile, LMatThres, PsiMatThres, 0);
%         FACluster_v0_3_stable_win(nFile+1, LMatThres, PsiMatThres, LifetimeThres);
%    end
% end

% % Visualization - Figure
Ablation_v_4_2(fishList, '_Double_FA_JointWin_DefinedSeg'); %plot change/ratio of synchronization level
% 
Ablation_v_5_1(fishList, '_definedSeg'); % plot experimental summary in matrix

Ablation_v_6(fishListAll, '_allFish'); % plot factor size vs. region length
