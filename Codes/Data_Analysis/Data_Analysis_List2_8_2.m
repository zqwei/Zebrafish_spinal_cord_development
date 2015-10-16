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


function Data_Analysis_List2_8_2(nFile)

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
    Thres               = 0.3; % This thres is used to uncover the final two groups
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
        
    figure;        
    h                   = suptitle(fileName);
    set(h,'Interpreter','none');  
        
    for nTime           = 1:length(timePoints)-1
        subplot(mRow, mCol, nTime)
        hold on;
        xTrack          = squeeze(mean(tracks(:, timePoints(nTime)+1:timePoints(nTime+1), 1), 2)); %#ok<NODEF>
        yTrack          = squeeze(mean(tracks(:, timePoints(nTime)+1:timePoints(nTime+1), 2), 2));
        LMatnTime       = LMatTime{nTime};
        hold on;        
        sizeCluster     = sum(abs(LMatTime{nTime})>Thres, 1);
        clusterToShow   = sizeCluster>minSizeCluster;
        LMatnTime       = abs(LMatnTime(:, clusterToShow));        
        NeuronsInNoFactor = (sum(abs(LMatnTime)>Thres, 2) == 0);
        plot(xTrack(NeuronsInNoFactor), yTrack(NeuronsInNoFactor), '.', 'MarkerFaceColor',[0.5 0.5 0.5], 'MarkerEdgeColor', [0.5 0.5 0.5],'Markersize',12);        
        for nFactor     = 1: size(LMatnTime, 2)
            plot(xTrack(LMatnTime(:,nFactor)>Thres), yTrack(LMatnTime(:,nFactor)>Thres), '.','Markersize',12);
        end
        hold off;
        ylim([000 400])
        xlim([0 1600])
        box off
%         legend('boxoff','location','northwest');
        xlabel('R-C location (um)')
        ylabel('M-L location (um)')
        title(['Time from: ' num2str(timePoints(nTime)/4/3600, '%.2f') ' to ' num2str(timePoints(nTime+1)/4/3600, '%.2f') ' hr'])
    end
    setPrint(mCol*8, mRow*6, [plotDir, 'FALMatClusterThre05_', fileName], 'pdf')    
end