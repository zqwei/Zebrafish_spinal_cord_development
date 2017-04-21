%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.  Covariance analysis: Factor analysis -- number of island using non
% LONO methods with selected neurons
% 
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function FACluster_v0_2_2(nFile)    

    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>  
    load([tempDatDir, 'LONOLoading_' fileName, '.mat'], 'CorrectedLMat') 
    numTime           = size(CorrectedLMat, 1);
    
    LONOM             = zeros(numTime, 1);

    for nTime         = 1:numTime    
        LMat          = CorrectedLMat{nTime};
        LMatNeuron    = LMat>0;
        numFactor     = size(LMatNeuron, 2);
        for nFactor   = 1:size(LMatNeuron, 2)
            for mFactor = nFactor+1:size(LMatNeuron, 2)
                if all(ismember(find(LMatNeuron(:,nFactor))', find(LMatNeuron(:,mFactor))')) || all(ismember(find(LMatNeuron(:,mFactor))', find(LMatNeuron(:,nFactor))'))
                    numFactor = numFactor - 1;
                    continue;
                end
            end
        end
        LONOM(nTime) = numFactor;
    end
    
%     LONOMLow          = zeros(numTime, 1);    
%     EVLONOMat         = squeeze(mean(EVLONO, 3));
%     for nTime         = 1:numTime    
%         if sum(~isnan(EVLONOMat(nTime, :))) > 0
%             maxEV        = nanmax(EVLONOMat(nTime, :));
%             if maxEV>0
%                 LONOMLow(nTime) = find(EVLONOMat(nTime, :)>0.85*maxEV, 1, 'first');
%             end
%         end        
%     end
%     
%     LONOMUp           = zeros(numTime, 1);    
%     EVLONOMat         = squeeze(mean(EVLONO, 3));
%     for nTime         = 1:numTime    
%         if sum(~isnan(EVLONOMat(nTime, :))) > 0
%             maxEV        = nanmax(EVLONOMat(nTime, :));
%             if maxEV>0
%                 LONOMUp(nTime) = find(EVLONOMat(nTime, :)>0.95*maxEV, 1, 'first');
%             end
%         end        
%     end
%     
    timePoints        = (1:numTime)';
    figure;  
    plot(timePoints/60, LONOM,'ko')
    ylabel('# factors')
    xlabel('Time (hour)')
    set(gca, 'xtick', 1:6)
    xlim([0 numTime/60])
    ylim([0 8])
    box off


    setPrint(8, 6, [plotDir, 'numFactorLONOActiveNeuronsDropOverlap_', fileName]);
    
    close all
        
end