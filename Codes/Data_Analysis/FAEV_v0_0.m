%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.  Explained variance: Factor analysis -- corrected loading matrix 
% 
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function FAEV_v0_0(nFile)            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>    
    mnx               = [];
    load([tempDatDir, fileName, '.mat'], 'dff', 'timePoints', 'side', 'mnx'); 
    load([tempDatDir, 'LONOLoading_' fileName, '.mat'], 'CorrectedLMat', 'PsiMat'); 
    
    numTime           = length(CorrectedLMat); %#ok<USENS>
    numNeuron         = size(dff, 1);
    EVMat             = nan(numNeuron, numTime);
    
    for  nTime        = 1:numTime
        lambda        = CorrectedLMat{nTime};
        lambda(isnan(lambda)) = 0;
        psi           = PsiMat{nTime};
        psi(sum(abs(lambda), 2)==0) = 1;
        slicedDFF     = dff(:,timePoints(nTime)+1:timePoints(nTime)+1200); 
        slicedDFF     = bsxfun(@minus, slicedDFF, mean(slicedDFF,2));
        slicedDFF     = bsxfun(@rdivide, slicedDFF, std(slicedDFF,[],2))';                 
        EVMat(:, nTime) = LONOFASingleUnitEV(slicedDFF, lambda, psi);
    end
    
    for nNeuron       = 1:numNeuron
        evNeuron      = smooth(EVMat(nNeuron, :), 41);
        EVMat(nNeuron, evNeuron > EVMat(nNeuron, :)'*5) = nan;
    end
    
    m                 = ceil(sqrt(numNeuron));
    
    figure;
    
    for nNeuron       = 1:numNeuron
        subplot(m, m, nNeuron)
        timeIndex = ~isnan(EVMat(nNeuron, :));
        if side(nNeuron) == 1
            plot(find(timeIndex)/60, EVMat(nNeuron, timeIndex), 'o', 'color', [0    0.4470    0.7410]);
        else
            plot(find(timeIndex)/60, EVMat(nNeuron, timeIndex), 'o', 'color', [0.6350    0.0780    0.1840]);
        end
        
        if ~isempty(mnx)
            if mnx(nNeuron)
                title('MNX+ neuron')
            else
                title('MNX- neuron')
            end
        end
        
        ylim([0 1])
        xlim([0 numTime/60])        
        box off
    end
    
    save([tempDatDir, 'EV_' fileName, '.mat'], 'EVMat')
    
    
end