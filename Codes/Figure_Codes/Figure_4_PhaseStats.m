%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Figure 4 Phase statistics
% calculate statistics of LR phase delay by the end of imaging
% For each fish, median delay of last 30 time windows, between largest ensembles on L & R
%
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%


function medDelays = Figure_4_PhaseStats()
datasets = [3, 4, 7, 12, 10, 13, 16]; 
% datasets = 15;
numWindows = 20; % # of windows to count from the end
medDelays = nan(numel(datasets), 1);
addpath('../Func');
setDir;
for i =1:numel(datasets)
    nFile = datasets(i);
    fileName          = fileNames{nFile};
    load([tempDatDir, 'FALONO_', fileName, '.mat'], 'dumpDuplicatedFactorLONOM');
    load([tempDatDir, 'EvoLoading_' fileName, '_v2.mat'], 'networkMat');
    numTime       = size(networkMat, 1); %#ok<*NODEF>
    delays = nan(numWindows, 1);
    for j = 1:numWindows
        nTime = numTime - j+1;
        networkTime    = networkMat{nTime, 1}; %#ok<*USENS>
        networkCorr    = networkMat{nTime, 2};
        numUnitsFactor = [networkTime.nUnit];
        sideFactor     = [networkTime.side];
        lFactor = find(sideFactor==1 & numUnitsFactor==max(numUnitsFactor(sideFactor==1)));
        rFactor = find(sideFactor==2 & numUnitsFactor==max(numUnitsFactor(sideFactor==2)));
        fid = sort([lFactor, rFactor]);
        delayCorr = networkCorr.corrMat(fid(1)-1, fid(2)-1);
        delayTime = abs(networkCorr.delayMat(fid(1)-1, fid(2)-1));
        if delayTime<10 && delayCorr>0.3
            delays(j) = delayTime;
        end
    end
    medDelays(i) = nanmedian(delays);
end

end
