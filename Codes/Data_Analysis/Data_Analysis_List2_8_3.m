%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  Covariance analysis: factor analysis
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.4 FA vs islands -- clustering
%     following 2.7
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function Data_Analysis_List2_8_3(nFile)

    if nargin<1
        nFile = 1;
    end
        
    addpath('../Func');
    setDir;    
    fileDirName         = fileDirNames{nFile}; %#ok<NASGU,USENS>
    fileName            = fileNames{nFile}; %#ok<USENS>
    load([tempDatDir, fileName, '.mat']);
    load([tempDatDir, fileName, '_FAEV.mat'], 'opt2Dim');
    load([tempDatDir, fileName, '_LMatTime.mat'], 'LMatTime')
    
    mCol                = 8;
    mRow                = ceil((length(timePoints)-1)/mCol);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define the final two groups -- moving left and right
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    Thres               = 0.5; % This thres is used to uncover the final two groups
    minSizeCluster      = 3;
    FinalGroup          = false(size(LMatTime{1},1), 2);                                %#ok<USENS>
    % Group 1
    [~, TimeG1]         = max(cell2mat(arrayfun(@(x) ...
                          max(sum(LMatTime{x}>Thres,1)), 1:length(LMatTime), ...
                          'UniformOutput', false))); 
    [~, FAIndexG1]      = max(sum(LMatTime{TimeG1}>Thres,1));
    FinalGroup(:,1)     = LMatTime{TimeG1}(:,FAIndexG1)>Thres;                      
    % Group 2
    [~, TimeG2]         = max(cell2mat(arrayfun(@(x) ...
                          max(sum(LMatTime{x}(~FinalGroup(:,1),:)>Thres,1)), ...
                          1:length(LMatTime), ...
                          'UniformOutput', false)));
    [~, FAIndexG2]      = max(sum(LMatTime{TimeG2}(~FinalGroup(:,1),:)>Thres,1));
    FinalGroup(:,2)     = LMatTime{TimeG1}(:,FAIndexG2)>Thres;  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Left and right groups defines the color of each group
    % Neurons belong to left or right groups has a saturation is equal to
    % zero; otherwise its saturation is equal to one
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %    
    
    G1G2Component       = [];
    NumNonClusNeurons   = zeros(length(timePoints)-1,1);
    
    for nTime           = 1:length(timePoints)-1
        LMatnTime       = LMatTime{nTime};       
        sizeCluster     = sum(abs(LMatTime{nTime})>Thres, 1);
        clusterToShow   = sizeCluster>minSizeCluster;
        LMatnTime       = abs(LMatnTime(:, clusterToShow));        
        NeuronsInNoFactor = (sum(abs(LMatnTime)>Thres, 2) == 0);
        NumNonClusNeurons(nTime) = sum(NeuronsInNoFactor);
        for nFactor     = 1: size(LMatnTime, 2)
            LMatTimeFactorPos = LMatnTime(:,nFactor)>Thres;
            nNeuronPosInG1    = sum(LMatTimeFactorPos & FinalGroup(:,1))/sum(LMatTimeFactorPos);
            nNeuronPosInG2    = sum(LMatTimeFactorPos & FinalGroup(:,2))/sum(LMatTimeFactorPos);
            G1G2Component     = [G1G2Component; nNeuronPosInG1, nNeuronPosInG2, timePoints(nTime)]; %#ok<AGROW>
        end
    end
    
    figure;
    plot(timePoints(1:nTime)/4/3600, size(LMatnTime,1)-NumNonClusNeurons,'o');
    xlabel('Time (hour)')
    ylabel('tot. # neurons')
    xlim([0 timePoints(end)/4/3600]); %#ok<COLND>
    setPrint(6, 4.5, [plotDir, 'FALMatTotalNumNeuron_', fileName], 'pdf')
    
    figure;
    hold on
    plot([0 1],[0 1],'--k','linewid',1.0)
    scatter(G1G2Component(:,1), G1G2Component(:,2),[], G1G2Component(:,3)/4/3600, 'o','fill')
    hold off
    colorbar
    xlabel('Ipsi. Comp. #1')
    ylabel('Ipsi. Comp. #2')
    axis([0 1 0 1]);
    setPrint(8, 6, [plotDir, 'FALMatClusterC1C2_', fileName], 'pdf')    
end