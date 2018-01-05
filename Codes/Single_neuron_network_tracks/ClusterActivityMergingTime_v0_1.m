%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LMatCorrection_v_0_1
% PreLMatTracker_v_0_1
% Spectrogram_v0_3
% Spectrogram_v0_4
%
% Based on the merging cluster identified by ClusterActivityMergingTime_v0_0
%
% -------------------------------------------------------------------------
% 
% Ziqiang Wei
% weiz@janelia.hhmi.org
%
%

function ClusterActivityMergingTime_v0_1(nFile)
    addpath('../Func');
    setDir;    
    fileName   = fileNames{nFile};   %#ok<*USENS>
    load([tempDatDir, fileName, '.mat'], 'dff', 'timePoints', 'new_x', 'new_y', 'new_z')
    load([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_1.mat'], 'new_activeNeuronMat', 'oscNeuronMat') 
    load([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_2.mat'], 'preLMat', 'preLMatIndex', 'preLMatTime') 
    numTime           = size(new_activeNeuronMat, 2);
    numNeuron         = size(new_activeNeuronMat, 1);
    neuronFactorIndex = false(size(new_activeNeuronMat));
    timeWin           = 10;

    numFactor     = max(preLMatIndex);
    corr_mat      = nan(numFactor, 4);
    
    for nFactor   = 1:numFactor
        timeInd   = preLMatTime(:, preLMatIndex == nFactor);
        if max(timeInd) - min(timeInd) > 5
            LMat  = preLMat(:, preLMatIndex == nFactor & preLMatTime == max(timeInd)); 
            timeLength = max(timeInd) - min(timeInd) - 5;
            timeLength = min(numTime - max(timeInd), timeLength);
            if max(timeInd) < numTime
                postLMatIndex         = preLMatIndex(:, preLMatTime == max(timeInd)+1);
                postLMat              = preLMat(:, preLMatTime == max(timeInd)+1);
                [max_value, similarity_index] = max(double(LMat') * double(postLMat));
                if isempty(max_value); continue; end
                globalLMat            = preLMat(:, preLMatIndex == postLMatIndex(similarity_index) & preLMatTime == max(timeInd)); 
                globalLMat            = sum(globalLMat, 2) > 0;
                if max_value == 0 || sum(globalLMat)==0 || sum(LMat)==0; continue; end
                localPreDff           = dff(LMat, timePoints(min(timeInd))+1:timePoints(min(timeInd)+timeLength));
                globaPreDff           = dff(globalLMat, timePoints(min(timeInd))+1:timePoints(min(timeInd)+timeLength));
                localPostDff          = dff(LMat, timePoints(max(timeInd))+1:timePoints(max(timeInd)+timeLength));
                globaPostDff          = dff(globalLMat, timePoints(max(timeInd))+1:timePoints(max(timeInd)+timeLength));
                corr_mat(nFactor, 1)  = sum(sum(abs(corr(localPreDff'))))/sum(LMat)/(sum(LMat)-1) - 1/(sum(LMat)-1);
                corr_mat(nFactor, 2)  = sum(sum(abs(corr(localPostDff'))))/sum(LMat)/(sum(LMat)-1) - 1/(sum(LMat)-1);
                corr_mat(nFactor, 3)  = sum(sum(abs(corr(globaPreDff'))))/sum(globalLMat)/(sum(globalLMat)-1) - 1/(sum(globalLMat)-1);
                corr_mat(nFactor, 4)  = sum(sum(abs(corr(globaPostDff'))))/sum(globalLMat)/(sum(globalLMat)-1) - 1/(sum(globalLMat)-1);
            end
        end
    end   
    figure;
    scatter(corr_mat(:,1), corr_mat(:,3), corr_mat(:,2)*200, corr_mat(:,4), 'filled')
    caxis([0 1])
    xlim([0 1])
    ylim([0 1])
    box off
    setPrint(8, 6, [plotNetDir 'LocalCommunityMergingCorr_' fileName], 'pdf')
    save([tempDatNetDir, 'LocalCommunityMerging_' fileName, '_v_0_1.mat'], 'corr_mat')
end