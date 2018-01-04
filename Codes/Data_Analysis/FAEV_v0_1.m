%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.  Explained variance: Factor analysis -- corrected loading matrix --
% sigmoid fit
% 
% computing half EV time
% computing half active time
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

function FAEV_v0_1(nFile)            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>    
    mnx               = [];
    load([tempDatDir, fileName, '.mat'], 'mnx', 'tracks', 'side', 'leafOrder', 'slicedIndex', 'activeNeuronMat'); 
    load([tempDatDir, 'EV_' fileName, '.mat'], 'EVMat')
    
    numTime           = size(EVMat, 2); %#ok<*NODEF>
    numNeuron         = size(EVMat, 1);
    halfEVTime        = nan(numNeuron, 1);
    RSquare           = zeros(numNeuron, 1);
    halfActTime       = nan(numNeuron, 1);
    firstActTime      = nan(numNeuron, 1);
    validFitIndex     = false(numNeuron, 1);
    
    m                 = ceil(sqrt(numNeuron));
    
    neuronName        = find(slicedIndex);
    neuronName        = neuronName(leafOrder);
    
    figure;
    timeBin           = 15;
    activeThres       = 0.65;
    
    for nNeuron       = 1:numNeuron
        subplot(m, m, nNeuron)
        
        if side(nNeuron) == 1
            mColor  = [0    0.4470    0.7410];
        else
            mColor  = [0.6350    0.0780    0.1840];
        end
        
        timeIndex = ~isnan(EVMat(nNeuron, :));
        actCurrNeuron   = activeNeuronMat(nNeuron, :);
        firstActiveTime = find(smooth(double(actCurrNeuron), timeBin) > activeThres, 1, 'first');
        if ~isempty(firstActiveTime)
            firstActTime(nNeuron) = firstActiveTime;
            removeIndex     = find(actCurrNeuron==0);
            removeIndex(removeIndex<firstActiveTime) = [];
            removeTimeIndex = false(1, numTime);
            removeTimeIndex(removeIndex) = true; % remove all zeros after first activation time
            timeIndex = timeIndex & ~ removeTimeIndex;
        else
            removeTimeIndex = false(1, numTime);
            timeIndex       = timeIndex & ~ removeTimeIndex;
        end
        
        if ~isempty(firstActiveTime)
            mdl = fitglm(find(~removeTimeIndex)/60, activeNeuronMat(nNeuron, ~removeTimeIndex), 'linear', 'Distribution', 'binomial');
            beta = mdl.Coefficients.Estimate;
            halfActTime(nNeuron) = -beta(1)/beta(2);
            if (mean(activeNeuronMat(nNeuron, 1:ceil(numTime*.5))) > 0.80 ||  mean(activeNeuronMat(nNeuron, :)) > 0.80) ...
                    && (halfActTime(nNeuron)<0 || halfActTime(nNeuron)>numTime/60)
                halfActTime(nNeuron) = 0;
            end
            if (halfActTime(nNeuron)<0 || halfActTime(nNeuron)>numTime/60); halfActTime(nNeuron) = nan; end
        end
                
        if sum(timeIndex) > 10 && halfActTime(nNeuron)>=0
            paddingLength            = 20;
            zerosPadding             = zeros(1, paddingLength+1);
            init_params              = [0, quantile(EVMat(nNeuron, timeIndex), 0.95), halfActTime(nNeuron), 1];
            
            if halfActTime(nNeuron) > numTime
                init_params(3)       = numTime;
                init_params(4)       = 0.001;
            end
            
            [fitParams, fitResult]   = sigm_fit([(-paddingLength:0)/60, find(timeIndex)/60], [zerosPadding, EVMat(nNeuron, timeIndex)], [0, nan, nan, nan], init_params, false);
            ypred                    = fitResult.ypred(paddingLength+2: end);
            ypredlowerCI             = fitResult.ypredlowerCI(paddingLength+2:end);
            ypredupperCI             = fitResult.ypredupperCI(paddingLength+2:end);
            RSquare(nNeuron)         = 1 - mean((EVMat(nNeuron, timeIndex)' - ypred).^2)./var(EVMat(nNeuron, timeIndex)'); %#ok<UDIM>
            halfEVTime(nNeuron)      = fitParams(3);
            if halfEVTime(nNeuron) < 0 || halfEVTime(nNeuron) > numTime/60
                halfEVTime(nNeuron)  = nan;
            end
%             [~, fitResult]           = sigm_fit(find(timeIndex)/60, EVMat(nNeuron, timeIndex)');
%             RSquare(nNeuron)         = 1 - mean((EVMat(nNeuron, timeIndex)' - fitResult.ypred).^2)./var(EVMat(nNeuron, timeIndex)');
        end
        

        
        hold on;
        plot(find(timeIndex)/60, EVMat(nNeuron, timeIndex), 'o', 'color', mColor);
%         plot((1:numTime)/60, activeNeuronMat(nNeuron, :), '.k')
        if ~isnan(halfEVTime(nNeuron)) && sum(activeNeuronMat(nNeuron, ~removeTimeIndex))>numTime*0.1
            validFitIndex(nNeuron) = true;
            plot(find(timeIndex)/60, ypred, '-', 'linewid', 2.0, 'Color', mColor);
            plot(find(timeIndex)/60, ypredlowerCI, '-', 'linewid', 0.5, 'Color', mColor); 
            plot(find(timeIndex)/60, ypredupperCI, '-', 'linewid', 0.5, 'Color', mColor); 
            % plot(find(~removeTimeIndex)/60, mdl.Fitted.Probability, '-k')
            
            title(['#' num2str(neuronName(nNeuron)) ' R^2=' num2str(RSquare(nNeuron))])
            if ~isempty(mnx)
                if mnx(nNeuron)
                    title({['#' num2str(neuronName(nNeuron)) ' MNX+; R^2=' num2str(RSquare(nNeuron), '%.2f')]; ['halfEV:' num2str(halfEVTime(nNeuron), '%.2f') 'halfAct:' num2str(halfActTime(nNeuron), '%.2f')]})
                else
                    title({['#' num2str(neuronName(nNeuron)) ' MNX-; R^2=' num2str(RSquare(nNeuron), '%.2f')]; ['halfEV:' num2str(halfEVTime(nNeuron), '%.2f') 'halfAct:' num2str(halfActTime(nNeuron), '%.2f')]})
                end
            end
        end
        
        hold off;
        ylim([0 1])
        xlim([0 numTime/60])    
        xlabel('Time (hour)')
        ylabel('EV')
        box off
    end
    
    setPrint(8*m, 6*m, [plotDir, 'SigFitEVActiveNeuron_', fileName], 'pdf');
    save([tempDatDir, 'EV_' fileName, '.mat'], 'halfActTime', 'halfEVTime', 'firstActTime', 'RSquare', 'validFitIndex', '-append');
end