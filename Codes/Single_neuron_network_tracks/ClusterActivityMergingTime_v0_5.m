%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LMatCorrection_v_0_1
% PreLMatTracker_v_0_1
% Spectrogram_v0_3
% Spectrogram_v0_4
%
% Based on the merging cluster identified by ClusterActivityMergingTime_v0_0
% Based on the merging cluster identified by
% ClusterActivityMergingTime_v0_1 with dense sampling every minutes
%
% -------------------------------------------------------------------------
% 
% Ziqiang Wei
% weiz@janelia.hhmi.org
%
%

addpath('../Func');
setDir;    

markerStyle = {'', '', 'o', '', '', '', 's', '', '', '^', '', '', '', '', '', '>'};

% % figure;
% % hold on
% % 
% % for nFile = [3 7 10 16]
% % 
% %     fileName   = fileNames{nFile};   %#ok<*USENS>
% %     load([tempDatDir, fileName, '.mat'], 'dff', 'timePoints', 'new_x', 'new_y', 'new_z')
% %     load([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_1.mat'], 'new_activeNeuronMat', 'oscNeuronMat') 
% %     load([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_2.mat'], 'preLMat', 'preLMatIndex', 'preLMatTime') 
% %     numTime           = size(new_activeNeuronMat, 2);
% %     numNeuron         = size(new_activeNeuronMat, 1);
% %     neuronFactorIndex = false(size(new_activeNeuronMat));
% %     numFactor     = max(preLMatIndex);
% %     
% %     for nFactor   = 1:numFactor
% %         timeInd   = preLMatTime(:, preLMatIndex == nFactor);        
% %         if max(timeInd) - min(timeInd) > 5 && max(timeInd) > 35 && max(timeInd) < numTime - 10
% %             merging_time = max(timeInd);
% %             % corr_mat is computed from merging time - 30 to merging time +
% %             % 30            
% %             LMat  = preLMat(:, preLMatIndex == nFactor & preLMatTime == merging_time); 
% %             postLMatIndex         = preLMatIndex(:, preLMatTime == merging_time+1);
% %             postLMat              = preLMat(:, preLMatTime == max(timeInd)+1);
% %             [max_value, similarity_index] = max(double(LMat') * double(postLMat));
% %             globalLMat            = preLMat(:, preLMatIndex == postLMatIndex(similarity_index) & preLMatTime == merging_time); 
% %             globalLMat            = sum(globalLMat, 2) > 0;
% %             if max_value == 0; continue; end
% %             corr_mat              = [];
% %             for nTime             = (merging_time - 10) : 5:  min(merging_time + 30, numTime)
% %                 localDff          = dff(LMat, timePoints(nTime)+(1:1200));
% %                 globalDff         = dff(globalLMat, timePoints(nTime)+(1:1200));
% %                 corr_local        = sum(sum(abs(corr(localDff'))))/sum(LMat)/(sum(LMat)-1) - 1/(sum(LMat)-1);
% %                 corr_global       = sum(sum(abs(corr(globalDff'))))/sum(globalLMat)/(sum(globalLMat)-1) - 1/(sum(globalLMat)-1);
% %                 corr_mat          = [corr_mat; corr_local, corr_global];
% %             end
% %         end
% %         
% % %         plot(corr_mat([1 3 end], 1), corr_mat([1 3 end], 2), 'linewid', 2, 'color', 'k');
% % %         plot(corr_mat(3, 1), corr_mat(3, 2), markerStyle{nFile}, 'linewid', 2, 'markerSize', 10);
% %         scatter(corr_mat([1 3 end], 1), corr_mat([1 3 end], 2), [], [1 2 3]', markerStyle{nFile});
% %     end   
% %         
% % end
% % 
% % xlim([0 1])
% % ylim([0 1])
% % plot([0 1], [0 1], '--k')
% % box off
% % set(gca, 'TickDir', 'out')
% % xlabel('Local correlation')
% % ylabel('Global correlation')
% % setPrint(8, 6, 'Real_signal_corr', 'pdf')
% % 


figure;
hold on

for nFile = [3 7 10 16]

    fileName   = fileNames{nFile};   %#ok<*USENS>
    load([tempDatDir, fileName, '.mat'], 'dff', 'timePoints', 'new_x', 'new_y', 'new_z')
    load([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_1.mat'], 'new_activeNeuronMat', 'oscNeuronMat') 
    load([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_2.mat'], 'preLMat', 'preLMatIndex', 'preLMatTime') 
    numTime           = size(new_activeNeuronMat, 2);
    numNeuron         = size(new_activeNeuronMat, 1);
    neuronFactorIndex = false(size(new_activeNeuronMat));
    numFactor         = max(preLMatIndex);
    
    dff_new           = nan(numNeuron, timePoints(end));
    
    for nNeuron   = 1:numNeuron
        dff_new(nNeuron, :) = ca2spike(dff(nNeuron, 1:timePoints(end)), new_activeNeuronMat(nNeuron, :));
    end
    
    for nFactor   = 1:numFactor
        timeInd   = preLMatTime(:, preLMatIndex == nFactor);        
        if max(timeInd) - min(timeInd) > 5 && max(timeInd) > 35 && max(timeInd) < numTime - 10
            merging_time = max(timeInd);
            % corr_mat is computed from merging time - 30 to merging time +
            % 30            
            LMat  = preLMat(:, preLMatIndex == nFactor & preLMatTime == merging_time); 
            postLMatIndex         = preLMatIndex(:, preLMatTime == merging_time+1);
            postLMat              = preLMat(:, preLMatTime == max(timeInd)+1);
            [max_value, similarity_index] = max(double(LMat') * double(postLMat));
            globalLMat            = preLMat(:, preLMatIndex == postLMatIndex(similarity_index) & preLMatTime == merging_time); 
            globalLMat            = sum(globalLMat, 2) > 0;
            if max_value == 0; continue; end
            corr_mat              = [];
            for nTime             = (merging_time - 10) : 5:  min(merging_time + 30, numTime)
                localDff          = dff_new(LMat, timePoints(nTime)+(1:1200));
                globalDff         = dff_new(globalLMat, timePoints(nTime)+(1:1200));
                corr_local        = abs(corr(localDff'));
                corr_local(eye(size(corr_local))==1) = nan;
                corr_local        = nanmean(corr_local(:));
                corr_global       = abs(corr(globalDff'));
                corr_global(eye(size(corr_global))==1) = nan;
                corr_global       = nanmean(corr_global(:));
                corr_mat          = [corr_mat; corr_local, corr_global];
            end
        end
        
%         plot(corr_mat([1 3 end], 1), corr_mat([1 3 end], 2), 'linewid', 2, 'color', 'k');
%         plot(corr_mat(3, 1), corr_mat(3, 2), markerStyle{nFile}, 'linewid', 2, 'markerSize', 10);
        corr_mat([1 3 end], :)
        scatter(corr_mat([1 3 end], 1), corr_mat([1 3 end], 2), [], [1 2 3]', markerStyle{nFile});
    end   
        
end

xlim([0 1])
ylim([0 1])
plot([0 1], [0 1], '--k')
box off
set(gca, 'TickDir', 'out')
xlabel('Local correlation')
ylabel('Global correlation')
setPrint(8, 6, 'Real_signal_corr', 'pdf')