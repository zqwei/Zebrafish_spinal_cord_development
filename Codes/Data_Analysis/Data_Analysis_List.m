% % function Data_Analysis_List_all
numFile = 3;

for nFile = 1:numFile        
%     disp(['Generate the orginal data file for dataset #' num2str(nFile)]);   
%     disp('==================================');
%     disp('Analysis #1 -- clustering analysis');   
%     disp('Generate clusters based on correlation matrix'); 
%     disp('Plot correlation matrix after being clustered');
%     disp('Check the colinearity after being clustered');
%     disp('Remove the time points with abrupt change of colinearity');
%     Neuron_selection_v0(nFile);

%     disp('==================================');
%     disp('Analysis #2 -- spectrogram');   
%     Spectrogram_v0_0(nFile) % Multitaper analysis
%     Spectrogram_v1_0(nFile) % Peak power cycle...
%     disp('Analysis #2 -- spectrogram -- test of white noise');   
%     Spectrogram_v1_1(nFile)

%     disp('==================================');
%     disp('Analysis #3 -- factor analysis -- number of factors');   
%     disp('Analysis #3 -- factor analysis -- number of factors -- non-LONO'); 
%     FACluster_v0_0(nFile) % # factors with all neurons
%     FACluster_v0_1(nFile) % # selected active neurons
%     FACluster_v0_2(nFile) % # factors with selected active neurons
%     disp('Analysis #3 -- factor analysis -- number of factors -- LONO'); 
%     FACluster_v1_0(nFile) % computing EV for n fold
%     FACluster_v1_1(nFile) % plot LONO - time-dim
%     FACluster_v1_2(nFile) % plot LONOM optimal number
%     FACluster_v2_0(nFile) % computing factors based on Varimax and side information
%     FACluster_v2_1(nFile) % thresholding and sorting, plot    
%     FACluster_v2_2(nFile) % missing neurons & missing neuron correction
%     FACluster_v2_3(nFile) % percentage pf active FA neurons vs # of factors
%     FACluster_v2_4(nFile) % number of neurons in each cluster
    
%     disp('==================================');
%     disp('--- Analysis # -- factor analysis -- factor contunity check');
%     disp('This is an older version of code, where all neurons are included in computation')
%     disp('Deprecate!')

%     disp('==================================');
%     disp('Analysis #4 -- factor analysis -- EV for single neurons'); %based on LONO
%     FAEV_v0_0(nFile)
%     FAEV_v0_1(nFile)


%     disp('==================================');
%     disp('Analysis #5 -- factor analysis -- Evolution of loading matrix');
%     disp('Analysis #5 -- factor analysis -- Evolution of loading matrix -- tree plot');
%     FACluster_v3_0(nFile)
%     FACluster_v3_1(nFile)

%     disp('==================================');
%     disp('Analysis #6 -- MNX');
%     MNX_v0_0(nFile)

%     disp('==================================');
%     disp('Analysis #7 -- Phase');
%     PhaseCluster_v0_0(nFile)


    disp('==================================');
    disp('Analysis #8 -- Analysis for all KS neurons');
%     FACluster_v0_3(nFile) % non-LONO number of factors
%     FACluster_v1_3(nFile) % EVLONO
%     FACluster_v1_4(nFile) % plot
%     FACluster_v1_5(nFile) % plot  
    FACluster_v2_5(nFile) % computing factors
    FACluster_v2_6(nFile) % correction factors
    FACluster_v3_2(nFile)
    MNX_v0_1(nFile)
    close all;
end

