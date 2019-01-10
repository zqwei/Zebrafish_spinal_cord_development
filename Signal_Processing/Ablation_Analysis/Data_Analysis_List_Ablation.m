% % function Data_Analysis_List_Ablation
addpath('../Data_Analysis/')

fishListCutA = [4, 3, 7, 2]; % anterior cut
fishListCutM = [1, 21, 9, 10, 12, 8]; % middle cut
fishListCutP = [16, 17, 13, 15]; % posterior cut
fishList = [fishListCutA, fishListCutM, fishListCutP];
fishList = 17;


for i = 1:numel(fishList)
    nFish = fishList(i);
    for nFile = [24+nFish*2-1];
    disp(['Generate the orginal data file for dataset #' num2str(nFile)]);
    if mod(nFile, 2) == 0 % after
        activeNeuronMat = Neuron_selection_v0_short_win(nFile-1);
        Neuron_selection_v0_1_short_win(nFile, activeNeuronMat);
    else % before
        activeNeuronMat = Neuron_selection_v0_short_win(nFile);
    end
    Neuron_selection_v1(nFile); % Correlation matrix using AP information
    Neuron_selection_v2(nFile); 
    FACluster_v0_0(nFile) % # factors with selected active neurons
    FACluster_v0_2_short_win(nFile) % plot LONOM optimal number -- LONO or non-LONO
    FACluster_v0_3_short_win(nFile) % compute loading matrix for optimal number FA
    FACluster_v0_3_whole_win(nFile)
    FACluster_v0_5_1_short_win(nFile) % generate FA evolution video
    close all;
    end
end

Ablation_v_2(fishListCutA, fishListCutM, fishListCutP);
Ablation_v_3_1(fishListCutA, fishListCutM, fishListCutP);