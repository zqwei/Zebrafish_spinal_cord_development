%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  Covariance analysis: factor analysis
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.4 FA vs islands -- lineage plot -- normalization
%     following 2.7
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function Data_Analysis_List2_9_1(nFile)

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
    
    plotScale           = 100;
    
%     numFactor           = zeros(length(timePoints)-1,1);
    
    figure;        
    title(fileName,'Interpreter','none');
    hold on;        
    for nTime           = 1:length(timePoints)-1
        LMatnTime       = LMatTime{nTime};
        LMatnTime       = LMatnTime(neuronIndex,:);
        NeuronsInNoFactor = (sum(abs(LMatnTime)>Thres, 2) == 0);
        drawCircle(find(NeuronsInNoFactor)*plotScale, nTime*plotScale*ones(sum(neuronIndex),1),...
            max(abs(LMatnTime(NeuronsInNoFactor,:)),[],2)*plotScale/2, [0 0 0]);
        
%         numFactor(nTime)= size(LMatnTime, 2);
        for nFactor     = 1:size(LMatnTime, 2)
            LMatTimeFactorPos = LMatnTime(:,nFactor)>Thres;
            nNeuronPosInG1    = sum(LMatTimeFactorPos & FinalGroup(:,1))/sum(LMatTimeFactorPos);
            nNeuronPosInG2    = sum(LMatTimeFactorPos & FinalGroup(:,2))/sum(LMatTimeFactorPos);
            LMatTimeFactorNeg = LMatnTime(:,nFactor)<-Thres;
            nNeuronNegInG1    = sum(LMatTimeFactorNeg & FinalGroup(:,1))/sum(LMatTimeFactorNeg);
            nNeuronNegInG2    = sum(LMatTimeFactorNeg & FinalGroup(:,2))/sum(LMatTimeFactorNeg);
            
            
            if nNeuronPosInG1+nNeuronPosInG2 > 0
                nNeuronPosInG1= nNeuronPosInG1/(nNeuronPosInG1+nNeuronPosInG2);
                nNeuronPosInG2= nNeuronPosInG2/(nNeuronPosInG1+nNeuronPosInG2);
                hValue        = (nNeuronPosInG1*0+nNeuronPosInG2*1)*0.5;
                sizeValue     = LMatnTime(LMatTimeFactorPos,nFactor);
                drawCircle(find(LMatTimeFactorPos)*plotScale, nTime*plotScale*ones(sum(LMatTimeFactorPos),1),...
                        sizeValue*plotScale/2, hsv2rgb(hValue, 1, 1));
            end
            
            if nNeuronNegInG1+nNeuronNegInG2 > 0
                nNeuronNegInG1= nNeuronNegInG1/(nNeuronNegInG1+nNeuronNegInG2);
                nNeuronNegInG2= nNeuronNegInG2/(nNeuronNegInG1+nNeuronNegInG2);                
                hValue        = 1 - (nNeuronNegInG1*0+nNeuronNegInG2*1)*0.5;
                sizeValue     = -LMatnTime(LMatTimeFactorNeg,nFactor);
                drawCircle(find(LMatTimeFactorNeg)*plotScale, nTime*plotScale*ones(sum(LMatTimeFactorNeg),1),...
                        sizeValue*plotScale/2, hsv2rgb(hValue, 1, 1));
            end
            
            if  nNeuronPosInG1+nNeuronPosInG2 == 0 && nNeuronNegInG1+nNeuronNegInG2 == 0
                sizeValue     = LMatnTime(LMatTimeFactorPos,nFactor);
                drawCircle(find(LMatTimeFactorPos)*plotScale, nTime*plotScale*ones(sum(LMatTimeFactorPos),1),...
                        sizeValue*plotScale/2, [0.2 0.2 0.2]); 
                sizeValue     = -LMatnTime(LMatTimeFactorNeg,nFactor);
                drawCircle(find(LMatTimeFactorNeg)*plotScale, nTime*plotScale*ones(sum(LMatTimeFactorNeg),1),...
                        sizeValue*plotScale/2, [0.6 0.6 0.6]);
            end                
        end 
        
        text((size(LMatnTime, 1)+1)*plotScale, nTime*plotScale, ['Factor #= ', num2str(size(LMatnTime, 2))]);
        
    end
    hold off;
    xlim([0 (numUnit+1)*plotScale])
    ylim([0 length(timePoints)*plotScale])
    xTickLabel                 = arrayfun(@(x) num2str(x,'%02d'), neuronIndex,'Uniform',false);
    yTickLabel                 = arrayfun(@(x) num2str(x,'%.2f'), timePoints(2:end)/4/3600,'Uniform',false);
    xlabel('Unit index');
    ylabel('Time (hour)');
    set(gca, 'XTick', (1:numUnit)*plotScale, 'XTickLabel', xTickLabel);
    set(gca, 'YTick', (1:length(timePoints)-1)*plotScale, 'YTickLabel', yTickLabel);
    
    setPrint(32, 24, [plotDir, 'FALineage_', fileName], 'pdf')
end