%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0.  Prepare statistics for leader cell analysis
%
% mnx_level_func
% factorSize: mean size of first active timeWindow=20 with outliers (>3sigma
% away) removed
% birthtime, siblins (if available)
% islet (if available)
%
% prepare a list of variable names that can be used for plotting
% 
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
% 

function Leader_v0(nFile) 
addpath('../Func');
setDir;
fileName          = fileNames{nFile};

% calculate factor size
windowSize = 20;
load([tempDatDir, fileName, '.mat'], 'timePoints', 'slicedIndex', 'leafOrder', 'activeNeuronMat');
load([tempDatDir, 'LONOLoading_', fileName, '.mat'], 'CorrectedLMat');
listLeaderMetrics = {'halfActTime', 'halfEVTime', 'diffHalfTime', 'factorSize', 'mnx', 'mnxFunc'};

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
    ft = factorSizeTimeSeries(i, ~isnan(factorSizeTimeSeries(i, :)));
    if ~isempty(ft) && numel(ft)>windowSize
        f = ft(1:windowSize);
        f(abs(f-mean(f))>3*std(f)) = [];
        factorSize(i) = mean(f);
    end
end
plotFactorSize(factorSizeTimeSeries, factorSize, plotDir, fileName);
save([tempDatDir, 'Leader_', fileName, '.mat'], 'factorSize');


% calculate support data
% mnx level at the end of functional imaging
load([fileDirNames{nFile} '\profile.mat'], 'mnx_level_func', 'birthtime', 'islet');
mnxFunc = mnx_level_func(slicedIndex, :);
mnxFunc = mnxFunc(leafOrder, 1);
save([tempDatDir, 'Leader_', fileName, '.mat'], 'mnxFunc', '-append');

% birthtime
if (exist('birthtime', 'var'))
    birthtime = birthtime(slicedIndex);
    birthtime = birthtime(leafOrder);
    save([tempDatDir, 'Leader_', fileName, '.mat'], 'birthtime', '-append');
    listLeaderMetrics{end+1} = 'birthtime';
end

if (exist('islet', 'var'))
    islet = islet(slicedIndex);
    islet = islet(leafOrder);
    save([tempDatDir, 'Leader_', fileName, '.mat'], 'islet', '-append');
    listLeaderMetrics{end+1} = 'islet';
end

save([tempDatDir, 'Leader_', fileName, '.mat'], 'listLeaderMetrics', '-append');
end

function plotFactorSize(factorSizeTimeSeries, factorSize, plotDir, fileName)
nNeurons = size(factorSizeTimeSeries, 1);
timePoints = (1:size(factorSizeTimeSeries, 2))*240;
nCol = 8;
nRow = ceil(nNeurons/nCol);
figure,

for i = 1:nNeurons
    subplot(nRow, nCol, i);
    hold on
    title(['#' num2str(i) ' factorSize:' num2str(factorSize(i), '%.2f')]);
    scatter(timePoints(factorSizeTimeSeries(i, :)>2)/3600/4, factorSizeTimeSeries(i, factorSizeTimeSeries(i, :)>2), 5, 'b');
    scatter(timePoints(factorSizeTimeSeries(i, :)==1)/3600/4, factorSizeTimeSeries(i, factorSizeTimeSeries(i, :)==1), 5, 'r');
    scatter(timePoints(factorSizeTimeSeries(i, :)==2)/3600/4, factorSizeTimeSeries(i, factorSizeTimeSeries(i, :)==2), 5, 'g');
    ylim([0 max(factorSizeTimeSeries(:))]);
    xlim([0 max(timePoints)/3500/4]);
    xlabel('Time (h)')
    ylabel('FactorSize')
    hold off
end

    setPrint(8*nCol, 6*nRow, [plotDir, 'FactorSizeEvolution_', fileName], 'pdf');

end