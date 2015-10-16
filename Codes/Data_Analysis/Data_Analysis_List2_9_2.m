%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  Covariance analysis: factor analysis
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.4 FA vs islands -- lineage plot -- bilateral index
%     following 2.9.1
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function Data_Analysis_List2_9_2(nFile)

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

    Thres               = 0.5; % This thres is used to uncover the final two groups
    numUnit             = size(LMatTime{1},1); %#ok<USENS>
    FinalGroup          = false(numUnit, 2);                                
    % Group 1
    [~, TimeG1]         = max(cell2mat(arrayfun(@(x) ...
                          max(sum(LMatTime{x}>Thres,1)), 1:length(LMatTime), ...
                          'UniformOutput', false))); 
    [~, FAIndexG1]      = max(sum(LMatTime{TimeG1}>Thres,1));
    FinalGroup(:,1)     = LMatTime{TimeG1}(:,FAIndexG1)>Thres;     
    numG1               = sum(FinalGroup(:,1));
    % Group 2
    [~, TimeG2]         = max(cell2mat(arrayfun(@(x) ...
                          max(sum(LMatTime{x}(~FinalGroup(:,1),:)>Thres,1)), ...
                          1:length(LMatTime), ...
                          'UniformOutput', false)));
    [~, FAIndexG2]      = max(sum(LMatTime{TimeG2}(~FinalGroup(:,1),:)>Thres,1));
    FinalGroup(:,2)     = LMatTime{TimeG1}(:,FAIndexG2)>Thres; 
    numG2               = sum(FinalGroup(:,2));
    
    neuronIndex         = 1:numUnit;
    neuronIndex         = [neuronIndex(FinalGroup(:,1)),neuronIndex(FinalGroup(:,2)),neuronIndex(~FinalGroup(:,1) & ~FinalGroup(:,2))];
    FinalGroup(:,1)     = [true(numG1,1);false(numUnit-numG1,1)];
    FinalGroup(:,2)     = [false(numG1,1);true(numG2,1);false(numUnit-numG1-numG2,1)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Left and right groups defines the color of each group
    % Neurons belong to left or right groups has a saturation is equal to
    % zero; otherwise its saturation is equal to one
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %    
    
    
    figure;        
%     title(fileName,'Interpreter','none');
    hold on;        
    plot([timePoints(1) timePoints(end)],[0.5 0.5],'--k','linewid',0.5) %#ok<COLND>
    for nTime           = 1:length(timePoints)-1
        LMatnTime       = LMatTime{nTime};
        LMatnTime       = LMatnTime(neuronIndex,:);
        totUnit         = size(LMatnTime,1);
        
        for nFactor     = 1:size(LMatnTime, 2)
            LMatTimeFactorPos = LMatnTime(:,nFactor)>Thres;
            nNeuronPosInG1    = sum(LMatTimeFactorPos & FinalGroup(:,1))/sum(LMatTimeFactorPos);
            nNeuronPosInG2    = sum(LMatTimeFactorPos & FinalGroup(:,2))/sum(LMatTimeFactorPos);
            LMatTimeFactorNeg = LMatnTime(:,nFactor)<-Thres;
            nNeuronNegInG1    = sum(LMatTimeFactorNeg & FinalGroup(:,1))/sum(LMatTimeFactorNeg);
            nNeuronNegInG2    = sum(LMatTimeFactorNeg & FinalGroup(:,2))/sum(LMatTimeFactorNeg);
            
            
            if nNeuronPosInG1+nNeuronPosInG2 > 0
                nNeuronPosInG1= nNeuronPosInG1/(nNeuronPosInG1+nNeuronPosInG2);
                plot(timePoints(nTime)/4/3600,nNeuronPosInG1,'o','Markersize',sum(LMatTimeFactorPos)/totUnit*20,'linewid',1)
            end
            
            if nNeuronNegInG1+nNeuronNegInG2 > 0
                nNeuronNegInG1= nNeuronNegInG1/(nNeuronNegInG1+nNeuronNegInG2);
                plot(timePoints(nTime)/4/3600,nNeuronNegInG1,'o','Markersize',sum(LMatTimeFactorNeg)/totUnit*20,'linewid',1)
            end               
        end
    end
    hold off;
    xlim([timePoints(1) timePoints(end)]/4/3600) %#ok<COLND>
    ylim([0 1])
    xlabel('Time (hour)');
    ylabel('Bilateral connectivity index');
    setPrint(8, 6, [plotDir, 'BilateralConIndex_', fileName], 'pdf')
end