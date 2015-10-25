%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.  Covariance analysis: Factor analysis -- Loading matrix using LONO 
% number of factors -- Varimax -- correction of loading matrix
% 
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function FACluster_v2_1(nFile)            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>    
    load([tempDatDir, fileName, '.mat'], 'sideSplitter', 'side'); 
    load([tempDatDir, fileName, '_numFactorNoLONOActiveNeuronsSimplified.mat'], 'LONOM')
    load([tempDatDir, fileName, '_LONOLoading.mat'], 'LMat', 'PsiMat');
    
    numPlot           = length(LMat); %#ok<USENS>
    CorrectedLMat     = LMat;
%     maxFactor         = max(LONOM);
    factorLoadThres   = 0.3;
    
    for nPlot         = 1:numPlot    
        if LONOM(nPlot)>0
%             mPlot         = mPlot + 1;
            LMatnTime     = LMat{nPlot};
            nFactor       = size(LMatnTime, 2);
%             LMatnTime(LMatnTime<min(nanmax(LMatnTime,[], 2))) = 0;
            LMatnTime(LMatnTime<factorLoadThres) = 0;

            numNeuronFactor = sum(LMatnTime>0, 1);
            sideFactor      = repmat(side, [1, nFactor]);
            sideFactor(~((LMatnTime)<=0)) = nan;
            sideNeuronFactor= nanmean(sideFactor, 1);
            [~, sortFactor] = sortrows([numNeuronFactor', sideNeuronFactor'], [2 -1]);
            LMatnTime       = LMatnTime(:, sortFactor);
            CorrectedLMat{nPlot} = LMatnTime;
        end
    end

    save([tempDatDir, fileName, '_LONOLoading.mat'], 'CorrectedLMat', '-append');
    
%     figure;
%     
%     m                 = ceil(sqrt(sum(LONOM>0)));
%     mPlot             = 0;
%     
%     for nPlot         = 1:numPlot    
%         if LONOM(nPlot)>0
%             mPlot         = mPlot + 1;
%             LMatnTime     = LMat{nPlot};
%             nFactor       = size(LMatnTime, 2);
%             LMatnTime(LMatnTime<min(nanmax(LMatnTime,[], 2))) = 0;
% 
%             numNeuronFactor = sum(LMatnTime>0, 1);
%             sideFactor      = repmat(side, [1, nFactor]);
%             sideFactor(~((LMatnTime)<=0)) = nan;
%             sideNeuronFactor= nanmean(sideFactor, 1);
%             [~, sortFactor] = sortrows([numNeuronFactor', sideNeuronFactor'], [2 -1]);
%             
%             LMatnTime       = LMatnTime(:, sortFactor);
% 
%             CorrectedLMat{nPlot} = LMatnTime;
%             subplot(m, m, mPlot)
%             hold on;
%             imagesc(LMatnTime, [0 1])
%             title(num2str(nPlot));
%             colormap(gray)
%             plot([1 maxFactor], [sideSplitter sideSplitter],'--r')
%             hold off;
%             axis off
%             xlim([0.5 maxFactor-0.5])
%         end
%     end
    
end