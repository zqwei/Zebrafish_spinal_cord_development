%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.  Covariance analysis: factor analysis
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.4 FA vs islands -- lineage plot -- bilateral index -- half life time
%     following 2.9.1
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function Data_Analysis_List2_9_3(nFile)

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
    
    LMatBilateralIndex  = LMatTime;
    LMatBilateralCon    = LMatTime;
    
    for nTime           = 1:length(timePoints)-1
        LMatnTime       = LMatTime{nTime};
        LMatnTime       = LMatnTime(neuronIndex,:);
        [totUnit, totFactor] = size(LMatnTime);
        LMatnTimeBI     = nan(totUnit, totFactor*2);
        LMatnTimeCN     = nan(totUnit, totFactor*2);
        
        for nFactor     = 1:size(LMatnTime, 2)
            LMatTimeFactorPos = LMatnTime(:,nFactor)>Thres;
            nNeuronPosInG1    = sum(LMatTimeFactorPos & FinalGroup(:,1))/sum(LMatTimeFactorPos);
            nNeuronPosInG2    = sum(LMatTimeFactorPos & FinalGroup(:,2))/sum(LMatTimeFactorPos);
            LMatTimeFactorNeg = LMatnTime(:,nFactor)<-Thres;
            nNeuronNegInG1    = sum(LMatTimeFactorNeg & FinalGroup(:,1))/sum(LMatTimeFactorNeg);
            nNeuronNegInG2    = sum(LMatTimeFactorNeg & FinalGroup(:,2))/sum(LMatTimeFactorNeg);
            
            
            if nNeuronPosInG1+nNeuronPosInG2 > 0
                nNeuronPosInG1= nNeuronPosInG1/(nNeuronPosInG1+nNeuronPosInG2);
                LMatnTimeBI(LMatTimeFactorPos, nFactor*2-1) = nNeuronPosInG1;
                LMatnTimeCN(LMatTimeFactorPos, nFactor*2-1) = LMatnTime(LMatTimeFactorPos,nFactor);
            end
            
            if nNeuronNegInG1+nNeuronNegInG2 > 0
                nNeuronNegInG1= nNeuronNegInG1/(nNeuronNegInG1+nNeuronNegInG2);
                LMatnTimeBI(LMatTimeFactorNeg, nFactor*2) = nNeuronNegInG1;
                LMatnTimeCN(LMatTimeFactorNeg, nFactor*2) = -LMatnTime(LMatTimeFactorNeg,nFactor);
            end               
        end
        
        LMatBilateralIndex(nTime)  = {LMatnTimeBI};
        LMatBilateralCon(nTime)    = {LMatnTimeCN};        
    end
    
    BilateralIndexNeuron = nan(totUnit, length(timePoints)-1);
    
    for nUnit           = 1:totUnit
        for nTime       = 1:length(timePoints)-1
            [~, nFactor] = nanmax(LMatBilateralCon{nTime}(nUnit, :));
            BilateralIndexNeuron(nUnit, nTime) = LMatBilateralIndex{nTime}(nUnit, nFactor);
        end
    end
    
    figure;   
    h = suptitle(fileName);
    set(h,'Interpreter','none');
    
    mCol             = 8;
    mRow             = ceil(totUnit/mCol);
    
    for nUnit        = 1:totUnit
        subplot(mRow, mCol, nUnit)
        if idx(nUnit) == 1
            mColor   = [ 0    0.4470    0.7410];
        elseif idx(nUnit) == 2
            mColor   = [0.6350    0.0780    0.1840];
        elseif idx(nUnit) == 3
            mColor   = [0 0 0];
        end
        plot(timePoints(1:end-1)/4/3600, BilateralIndexNeuron(nUnit,:),'o','color',[0.5 0.5 0.5]); %#ok<COLND>
        fitIndex           = ~isnan(BilateralIndexNeuron(nUnit,:));
        if sum(fitIndex) >= 6
            [~, fitResult]     = sigm_fit((timePoints(fitIndex)/4/3600), BilateralIndexNeuron(nUnit,fitIndex));
            hold on;
            plot((timePoints(fitIndex)/4/3600), fitResult.ypred, '-', 'linewid', 2.0, 'Color', mColor);
            plot((timePoints(fitIndex)/4/3600), fitResult.ypredlowerCI, '-', 'linewid', 0.5, 'Color', mColor);
            plot((timePoints(fitIndex)/4/3600), fitResult.ypredupperCI, '-', 'linewid', 0.5, 'Color', mColor);
            hold off;
            RSquare             = 1 - mean((BilateralIndexNeuron(nUnit,fitIndex)' - fitResult.ypred).^2)./var(BilateralIndexNeuron(nUnit,fitIndex));        

        else
            RSquare             = 0;
        end
        ylim([0 1])
        xlim([timePoints(1)/4/3600 timePoints(end)/4/3600]) %#ok<COLND>
        title(['R^2 = ', num2str(RSquare,'%.3f')])    
        xlabel('Time (hour)')
        ylabel('Bilateral index')
        box off
    end
    setPrint(mCol*6, mRow*4.5, [plotDir, 'BilateralIndexNeuron_', fileName], 'pdf');
end