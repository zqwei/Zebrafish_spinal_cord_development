% % function Data_Analysis_List_Ablation
addpath('../Data_Analysis/')
% fishList = [19,22,11,25,26,20,24,23];


fishList = 30;


%% FA 
for i = 1:numel(fishList)
    for nFile = 23+fishList(i)*2 + [0, 1];
        disp(['Generate the orginal data file for dataset #' num2str(nFile)]);
        Neuron_selection_v0_short_win(nFile);

        Neuron_selection_v1(nFile); % Correlation matrix using AP information
        Neuron_selection_v2(nFile); 

        FACluster_v0_0(nFile) % # factors with selected active neurons
        FACluster_v0_2_short_win(nFile) % plot LONOM optimal number -- LONO or non-LONO
% 
        FACluster_v0_3_short_win(nFile) % compute loading matrix for optimal number FA
        FACluster_v0_4(nFile) % correct FA, removing the projection without low EV

    close all;
    end
end

% Analysis of FA result
LMatThres = 0.3;
PsiMatThres = 0.5;
LifetimeThres = 1; % factor life threshold only applies to after ablation
for i = 1:numel(fishList)
   for  nFile = 23+fishList(i)*2 +  1
        FACluster_v0_3_whole_win(nFile, LMatThres, PsiMatThres)
        FACluster_v0_5_1_short_win(nFile, LMatThres, PsiMatThres) % generate FA evolution video
        FACluster_v0_3_stable_win(nFile, LMatThres, PsiMatThres, 0);
        FACluster_v0_3_stable_win(nFile+1, LMatThres, PsiMatThres, LifetimeThres);
   end
end
