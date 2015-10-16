% % function Data_Analysis_List_all
numFile = 3;

for nFile = numFile %1:numFile        
    disp(['Generate the orginal data file for dataset #' num2str(nFile)]);   
    disp('Remove the units with strong correlation and close locations');  
    Data_Analysis_List0_1(nFile);
    
    disp('==================================');
    disp('Analysis #1 -- clustering analysis');   
    disp('Generate clusters based on correlation matrix'); 
    disp('Remove the units outside clusters with minimal size = 3');
%     Data_Analysis_List1(nFile);
    disp('Plot correlation matrix after being clustered');
    Data_Analysis_List1_1(nFile);
    disp('Check the colinearity after being clustered');
    disp('Remove the time points with abrupt change of colinearity');
    Data_Analysis_List1_2(nFile);
    
    disp('==================================');
    disp('Analysis #2 -- factor analysis');   
    disp('Compute EV as a funtion of time and the number of factor');
    Data_Analysis_List2(nFile);
    disp('Plot EV as a funtion of time and the number of factor');
    disp('Find optimal number of factor at each time point');
    disp('Check EV for each neuron');
    Data_Analysis_List2_1(nFile);
    
    disp('==================================');
    disp('Analysis #3 -- factor analysis -- contunity check');   
    disp('Compute EV as a funtion across time points');
    Data_Analysis_List2_2(nFile);
    disp('Plot EV as a funtion across time points');
    Data_Analysis_List2_3(nFile);
    
    disp('==================================');
    disp('Analysis #4 -- factor analysis -- summary');  
    Data_Analysis_List2_0(nFile); %FACrossTime
    disp('Analysis #4 -- factor analysis -- number of factors using multiple tests'); 
    Data_Analysis_List2_0_4(nFile); %Num of Cluster Cross Time Non LONO methods
    Data_Analysis_List2_0_3(nFile); %FANumClusterCrossTime
    
    disp('==================================');
    disp('Analysis #5 -- factor analysis -- EV for single neurons');
    % group is predefined by Yinan using PCA or FA according to 
    % last time point
    Data_Analysis_List2_0_1(nFile); % FAEVSingleUnit
    Data_Analysis_List2_0_2(nFile);
     % half time
    Data_Analysis_List2_8_4(nFile); % HalfTimeThre05
    Data_Analysis_List2_4(nFile);
    
    disp('==================================');
    disp('Analysis #6 -- factor analysis -- Computing loading matrix');
    Data_Analysis_List2_7_0(nFile);
    disp('Analysis #6 -- factor analysis -- Matrix decomposition');
    Data_Analysis_List2_7_1(nFile); % FADecompositionCrossTime
    Data_Analysis_List2_7_3(nFile); % FADecompositionPsiCrossTime
    disp('Analysis #6 -- factor analysis -- Matrix decomposition -- Rotation');
    Data_Analysis_List2_7_2(nFile); % FADecompositionRotation
    disp('Analysis #6 -- factor analysis -- Evolution of loading matrix');
%     Data_Analysis_List2_7(nFile);
    Data_Analysis_List2_5(nFile); % FALMat
    disp('Analysis #6 -- factor analysis -- Evolution of loading matrix -- number of neurons in each cluster');
    Data_Analysis_List2_5_1(nFile); 
    disp('Analysis #7 -- factor analysis -- Evolution of loading matrix -- Tree structure');
%     Data_Analysis_List2_6(nFile);
    disp('Analysis #8 -- factor analysis -- Evolution of groups');
    Data_Analysis_List2_8(nFile);
%     Data_Analysis_List2_8_1(nFile);
    Data_Analysis_List2_8_2(nFile);
    Data_Analysis_List2_8_3(nFile)
    disp('Analysis #9 -- factor analysis -- Lineage of cell');
    Data_Analysis_List2_9(nFile);
    Data_Analysis_List2_9_1(nFile);
    Data_Analysis_List2_9_2(nFile);
    Data_Analysis_List2_9_3(nFile);
    Data_Analysis_List2_9_4(nFile);
    Data_Analysis_List2_9_5(nFile);
    Data_Analysis_List2_9_6(nFile);
    
%     disp('==================================');
%     disp('Analysis #10 -- crossing-middle-line analysis -- Covariance matrix ordered by first few time points');
%     Data_Analysis_List3_1(nFile);
%     Data_Analysis_List3_2(nFile);
%     disp('Analysis #11 -- crossing-middle-line analysis -- Shuffling test of correlation matrix');
%     Data_Analysis_List3_3(nFile);

    close all;
end

