% % function Data_Analysis_List_all
numFile = 17;

for nFile = 1:numFile              
    disp(['Generate the orginal data file for dataset #' num2str(nFile)]);   
    disp('==================================');
    disp('Analysis #1 -- clustering analysis');   
    disp('Generate clusters based on correlation matrix'); 
    disp('Plot correlation matrix after being clustered');
    disp('Check the colinearity after being clustered');
    disp('Remove the time points with abrupt change of colinearity');
    Neuron_selection_v0(nFile);
    Neuron_selection_v1(nFile); % Correlation matrix using AP information    

    disp('==================================');
    disp('Analysis #2 -- factor analysis -- number of factors');   
    disp('Analysis #2 -- factor analysis -- number of factors -- non-LONO'); 
    disp('Analysis #2 -- factor analysis -- number of factors -- LONO'); 
    FACluster_v0_0(nFile) % # factors with selected active neurons
    FACluster_v0_1(nFile) % plot LONO - time-dim (pcolor & contourf)
    FACluster_v0_2(nFile) % plot LONOM optimal number -- LONO or non-LONO
    FACluster_v0_2_1(nFile) % plot LONOM optimal number with Confidence Interval
    FACluster_v0_3(nFile) % compute loading matrix for optimal number FA
    FACluster_v0_4(nFile) % correct FA, removing the projection without low EV
    FACluster_v0_5_1(nFile) % generate FA evolution video

    FAEV_v0_0(nFile) % generated EV time mat file
    FAEV_v0_1(nFile) % compute half EV time
    FAEV_v0_2(nFile) % compare half EV time vs activation time
    FAEV_v0_3(nFile) % EV time with location
    
    close all;
end

