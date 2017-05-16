%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LMatCorrection_v_0_1
% PreLMatTracker_v_0_1
% Spectrogram_v0_3
% Spectrogram_v0_4
%
% Track the endtime of local clusters
%
% -------------------------------------------------------------------------
% 
% Ziqiang Wei
% weiz@janelia.hhmi.org
%
%

function ClusterActivityMergingTime_v0_0(nFile)
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

    numFactor     = max(preLMatIndex);
    
    for nFactor   = 1:numFactor
        timeInd   = preLMatTime(:, preLMatIndex == nFactor);

        if max(timeInd) - min(timeInd) > 5
            LMat      = preLMat(:, preLMatIndex == nFactor & preLMatTime == max(timeInd));            
            if max(timeInd) < numTime
                postLMatIndex = preLMatIndex(:, preLMatTime == max(timeInd)+1);
                postLMat  = preLMat(:, preLMatTime == max(timeInd)+1);
                [max_value, similarity_index] = max(double(LMat') * double(postLMat));
                globalLMat = preLMat(:, preLMatIndex == postLMatIndex(similarity_index) & preLMatTime == max(timeInd)); 
                if max_value > 0
                    disp([nFile, nFactor, min(timeInd), max(timeInd), postLMatIndex(similarity_index), sum(LMat), sum(sum(globalLMat,2)>0)])
                end
            end
        end
    end   
        
end