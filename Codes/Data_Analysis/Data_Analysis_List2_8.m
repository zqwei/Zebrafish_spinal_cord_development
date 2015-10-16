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


function Data_Analysis_List2_8(nFile)

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
%     maxWeight            = cell2mat(arrayfun(@(x) ...
%                           max(LMatTime{x},[],2), 1:length(LMatTime), ...
%                           'UniformOutput', false));
%     minWeight            = cell2mat(arrayfun(@(x) ...
%                           min(LMatTime{x},[],2), 1:length(LMatTime), ...
%                           'UniformOutput', false));
%     imagesc(abs(minWeight)>abs(maxWeight))

    mCol                = 8;
    mRow                = ceil((length(timePoints)-1)/mCol);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define the final two groups -- moving left and right
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    Thres               = 0.3; % This thres is used to uncover the final two groups
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
    
    plotScale           = 100;
    
    figure;        
    h                   = suptitle(fileName);
    set(h,'Interpreter','none');  
        
    for nTime           = 1:length(timePoints)-1
        subplot(mRow, mCol, nTime)
        xTrack          = squeeze(mean(tracks(:, timePoints(nTime)+1:timePoints(nTime+1), 1), 2)); %#ok<NODEF>
        yTrack          = squeeze(mean(tracks(:, timePoints(nTime)+1:timePoints(nTime+1), 2), 2));
        LMatnTime       = LMatTime{nTime};
        hold on;
        NeuronsInNoFactor = (sum(abs(LMatnTime)>Thres, 2) == 0);
        scatter(xTrack(NeuronsInNoFactor), yTrack(NeuronsInNoFactor), max(abs(LMatnTime(NeuronsInNoFactor,:)),[],2)*plotScale,'MarkerFaceColor','k', 'MarkerEdgeColor', 'k');
        
        for nFactor     = size(LMatnTime, 2):-1:1
            LMatTimeFactorPos = LMatnTime(:,nFactor)>Thres;
            nNeuronPosInG1    = sum(LMatTimeFactorPos & FinalGroup(:,1))/sum(LMatTimeFactorPos);
            nNeuronPosInG2    = sum(LMatTimeFactorPos & FinalGroup(:,2))/sum(LMatTimeFactorPos);
            LMatTimeFactorNeg = LMatnTime(:,nFactor)<-Thres;
            nNeuronNegInG1    = sum(LMatTimeFactorNeg & FinalGroup(:,1))/sum(LMatTimeFactorNeg);
            nNeuronNegInG2    = sum(LMatTimeFactorNeg & FinalGroup(:,2))/sum(LMatTimeFactorNeg);
            
            if nNeuronPosInG1+nNeuronPosInG2 > 0
                hValue        = (nNeuronPosInG1*0+nNeuronPosInG2*1)*0.5;
                sizeValue     = LMatnTime(LMatTimeFactorPos,nFactor);
                scatter(xTrack(LMatTimeFactorPos), yTrack(LMatTimeFactorPos), sizeValue*plotScale,'MarkerFaceColor', hsv2rgb(hValue, 1, 1), 'MarkerEdgeColor', hsv2rgb(hValue, 1, 1));
            end
            
            if nNeuronNegInG1+nNeuronNegInG2 > 0
                hValue        = 1 - (nNeuronNegInG1*0+nNeuronNegInG2*1)*0.5;
                sizeValue     = -LMatnTime(LMatTimeFactorNeg,nFactor);
                scatter(xTrack(LMatTimeFactorNeg), yTrack(LMatTimeFactorNeg), sizeValue*plotScale,'MarkerFaceColor', hsv2rgb(hValue, 1, 1), 'MarkerEdgeColor', hsv2rgb(hValue, 1, 1)); 
            end
            
            if  nNeuronPosInG1+nNeuronPosInG2 == 0 && nNeuronNegInG1+nNeuronNegInG2 == 0
                sizeValue     = LMatnTime(LMatTimeFactorPos,nFactor);
                scatter(xTrack(LMatTimeFactorPos), yTrack(LMatTimeFactorPos), sizeValue*plotScale,'MarkerFaceColor', [0.3 0.3 0.3], 'MarkerEdgeColor', [0.3 0.3 0.3]); 
                sizeValue     = -LMatnTime(LMatTimeFactorNeg,nFactor);
                scatter(xTrack(LMatTimeFactorNeg), yTrack(LMatTimeFactorNeg), sizeValue*plotScale,'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', [0.7 0.7 0.7]); 
            end                
        end
        hold off;
        ylim([000 400])
        xlim([0 1600])
        box off
        title(['Time from: ' num2str(timePoints(nTime)/4/3600, '%.2f') ' to ' num2str(timePoints(nTime+1)/4/3600, '%.2f') ' hr'])
    end
    setPrint(mCol*16, mRow*12, [plotDir, 'FALMatCluster_', fileName], 'pdf')    
end