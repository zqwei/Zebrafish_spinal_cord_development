%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0.  Prepare statistics for leader cell analysis
%
% mnx_level_func
% factorSize: mean size of first active timeWindow=20 after firstActTime
% birthtime, siblins (if available)
% islet (if available)
%
% prepare a list of variable names that can be used for plotting
% 
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
% 

function [numOutlier, numFit] = Leader_v0(nFile, RSquareThres) 
addpath('../Func');
setDir;
fileName          = fileNames{nFile};


load([tempDatDir, fileName, '.mat'], 'timePoints', 'slicedIndex', 'leafOrder', 'activeNeuronMat');
load([tempDatDir, 'LONOLoading_', fileName, '.mat'], 'CorrectedLMat');
load([tempDatDir, 'EV_', fileName, '.mat'], 'halfActTime', 'halfEVTime', 'firstActTime', 'validFitIndex', 'RSquare', 'halfEVTimeCI');
listLeaderMetrics = {'halfActTime', 'halfEVTime', 'diffHalfTime', 'factorSize', 'mnx', 'mnxFunc', 'activeTime', 'patternTime'};

% pattern time and active time
patternTime = halfEVTime;
patternTime(RSquare<=RSquareThres | ~validFitIndex)= NaN;
patternTimeBound = halfEVTimeCI(:, 2);
patternTimeBound(RSquare<=RSquareThres | ~validFitIndex) = NaN;
activeTime = firstActTime/60;
numOutlier = sum((activeTime - patternTimeBound)>0);
numFit = sum(~isnan(patternTime) & halfEVTimeCI(:, 2)-halfEVTimeCI(:, 1)<1);
% save([tempDatDir, 'Leader_', fileName, '.mat'], 'patternTime', 'activeTime');

% calculate factor size
windowSize = 20;
nNeurons = numel(leafOrder);
factorSizeTimeSeries = nan(nNeurons, numel(timePoints));

for period = 1:numel(timePoints)
    LMat = CorrectedLMat{period};
    for i = 1:nNeurons
        if activeNeuronMat(i, period)
            % determine factor belongings - a single active neuron
            if any(isnan(LMat(i, :))) || ~any(LMat(i, :)>0)
                factorSizeTimeSeries(i, period) = 1;
            else
                % factor belongings - determine dominate factor
                factorSizeTimeSeries(i, period) = sum(LMat(:, LMat(i, :)==max(LMat(i, :)))>0);
            end
        end
    end
end

factorSize = nan(nNeurons, 1);
for i = 1:nNeurons
    if ~isnan(firstActTime(i))
        factorSizeTimeSeries(i, 1:firstActTime(i)-1) = NaN;
        ft = factorSizeTimeSeries(i, ~isnan(factorSizeTimeSeries(i, :)));
        if ~isempty(ft) && numel(ft)>windowSize
            f = ft(1:windowSize);
    %         f(abs(f-mean(f))>3*std(f)) = [];
            factorSize(i) = mean(f);
        end
    end
end


% % plot factor size over time
% nCol = 8;
% nRow = ceil(nNeurons/nCol);
% figure,
% 
% for i = 1:nNeurons
%     subplot(nRow, nCol, i);
%     hold on
%     titleText = sprintf(['#' num2str(i) ' FSize:' num2str(factorSize(i), '%.2f') '\nactTime:' num2str(activeTime(i), '%.2f')  ' patternTime:' num2str(patternTime(i), '%.2f')]);
%     title(titleText);
%     scatter(timePoints(factorSizeTimeSeries(i, :)>2)/3600/4, factorSizeTimeSeries(i, factorSizeTimeSeries(i, :)>2), 5, 'b');
%     scatter(timePoints(factorSizeTimeSeries(i, :)==1)/3600/4, factorSizeTimeSeries(i, factorSizeTimeSeries(i, :)==1), 5, 'r');
%     scatter(timePoints(factorSizeTimeSeries(i, :)==2)/3600/4, factorSizeTimeSeries(i, factorSizeTimeSeries(i, :)==2), 5, 'g');
% %     plot([halfActTime(i), halfActTime(i)], [0, max(factorSizeTimeSeries(:))], 'k-');
%     plot([activeTime(i), activeTime(i)], [0, max(factorSizeTimeSeries(:))], 'r-');
%     plot([patternTime(i), patternTime(i)], [0, max(factorSizeTimeSeries(:))], 'g-');
% %     plot([halfEVTime(i), halfEVTime(i)], [0, max(factorSizeTimeSeries(:))], 'g-');
%     ylim([0 max(factorSizeTimeSeries(:))]);
%     xlim([0 max(timePoints)/3600/4]);
%     xlabel('Time (h)')
%     ylabel('FactorSize')
%     hold off
% end
% setPrint(8*nCol, 6*nRow, [plotDir, 'FactorSizeEvolution_', fileName '_Rthres', num2str(RSquareThres)], 'pdf');
% close
%     
% save([tempDatDir, 'Leader_', fileName, '.mat'], 'factorSize', '-append');
% 
% 
% 
% % calculate support data
% % mnx level at the end of functional imaging
% load([fileDirNames{nFile} '\profile.mat'], 'mnx_level_func', 'birthtime', 'islet');
% if (exist('mnx_level_func', 'var'))
%     mnxFunc = mnx_level_func(slicedIndex, :);
%     mnxFunc = mnxFunc(leafOrder, 1);
%     save([tempDatDir, 'Leader_', fileName, '.mat'], 'mnxFunc', '-append');
% end
% 
% % birthtime
% if (exist('birthtime', 'var'))
%     birthtime = birthtime(slicedIndex);
%     birthtime = birthtime(leafOrder);
%     save([tempDatDir, 'Leader_', fileName, '.mat'], 'birthtime', '-append');
%     listLeaderMetrics{end+1} = 'birthtime';
% end
% 
% if (exist('islet', 'var'))
%     islet = islet(slicedIndex);
%     islet = islet(leafOrder);
%     save([tempDatDir, 'Leader_', fileName, '.mat'], 'islet', '-append');
%     listLeaderMetrics{end+1} = 'islet';
% end
% 
% save([tempDatDir, 'Leader_', fileName, '.mat'], 'listLeaderMetrics', '-append');
end
