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

function ClusterActivityMergingTime_v0_2(nFile)
    addpath('../Func');
    setDir;    
    fileName   = fileNames{nFile};   %#ok<*USENS>
    load([tempDatDir, fileName, '.mat'], 'dff', 'timePoints', 'new_x', 'new_y', 'new_z')
    load([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_1.mat'], 'new_activeNeuronMat', 'oscNeuronMat') 
    load([tempDatNetDir, 'LONOLoading_' fileName, '_v_0_2.mat'], 'preLMat', 'preLMatIndex', 'preLMatTime') 
    numTime           = size(new_activeNeuronMat, 2);
    numNeuron         = size(new_activeNeuronMat, 1);
    neuronFactorIndex = false(size(new_activeNeuronMat));
    timeWin           = 40;
    
    dffs          = dff;
    for nNeuron   = 1:numNeuron
        dffs(nNeuron, :) = ca2spike(dff(nNeuron, :), new_activeNeuronMat(nNeuron, :), timePoints);
    end

    numFactor     = max(preLMatIndex);
    corr_mat      = nan(numFactor, 3);
    
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
                globalLMat            = preLMat(:, preLMatIndex == postLMatIndex(similarity_index) & preLMatTime == max(timeInd)); 
                globalLMat            = sum(globalLMat, 2) > 0;
                if max_value == 0; continue; end
                localPreDff           = dffs(LMat, timePoints(min(timeInd)):timePoints(min(timeInd)+timeLength));
                localPreDff           = sum(localPreDff)>0;
                globaPreDff           = dffs(globalLMat, timePoints(min(timeInd)):timePoints(min(timeInd)+timeLength));
                globaPreDff           = sum(globaPreDff)>0;
                allPreDff             = (globaPreDff + localPreDff) >0;
                localPostDff          = dffs(LMat, timePoints(max(timeInd)):timePoints(max(timeInd)+timeLength));
                localPostDff          = sum(localPostDff)>0;
                globaPostDff          = dffs(globalLMat, timePoints(max(timeInd)):timePoints(max(timeInd)+timeLength));
                globaPostDff          = sum(globaPostDff)>0;
                allPostDff            = (globaPostDff + localPostDff) >0;
                corr_mat(nFactor, 1)  = max(xcorr(localPreDff', localPostDff', 'coeff'));
                corr_mat(nFactor, 2)  = max(xcorr(globaPreDff', globaPostDff', 'coeff'));
                corr_mat(nFactor, 3)  = max(xcorr(allPreDff', allPostDff', 'coeff'));
            end
        end
    end   
    figure;
    scatter(corr_mat(:,1), corr_mat(:,2), [], corr_mat(:,3), 'filled')
    caxis([0 1])
    xlim([0 1])
    ylim([0 1])
    box off
    setPrint(8, 6, [plotNetDir 'LocalCommunityMergingPrePostCorr_' fileName], 'pdf')
    corr_prepost_mat = corr_mat;
    save([tempDatNetDir, 'LocalCommunityMerging_' fileName, '_v_0_1.mat'], 'corr_prepost_mat', '-append')
end


function dffs_marks = ca2spike(dffs, activeNeuronTime, timePoints)
    dffs_mark = ones(240, 1) * activeNeuronTime;
    dffs_mark = dffs_mark(:);
    [peaks, locs] = findpeaks(dffs(1:timePoints(numel(timePoints))), 'minPeakHeight', 0.0, 'minPeakDistance', 8);
    dffs_mark_locs = dffs_mark(locs);
    locs = locs(peaks > max([peaks(dffs_mark_locs==0), 0]) + 0.01);
    dffs_marks = false(size(dffs));
    dffs_marks(locs) = true;
    dffs_marks = double(dffs_marks);
end