%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Signal Extraction Step 2: Pre-select neurons for curation
% Select neurons that meet at least one of the following criterion
% 1) positively-biased signal from Kolmogorov–Smirnov test 
% 2) Factor analysis result of timeRange with factor loading>FA_thres
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%
function Signal_Extraction_s2(inputFolder, FA_range, FA_thres)
addpath('../Func');
% adaptation: use a sliding window of 5 min data every 1 min. Keep all the neurons that is ever detected "active"
load([inputFolder '.simpleTrack\active_neuron_extraction\profile.mat']);
windowSize = 1200; % 5 min
timeInterval = 240; % 1min
backgroundLevel = 90;
active_thres = 4; % must be active for this number of time windows to qualify

% baseline estimation using sgolay filter
numOrder      = 9;
lenWindow     = 511;
baseline      = sgolayfilt(profile_all, numOrder, lenWindow, [], 2);

    
dff = bsxfun(@rdivide,(profile_all - baseline), nanmean(baseline, 2)-backgroundLevel);



activeNeuronMat = false(nCells, numel(timepoints)/timeInterval);
for t_start = 1:timeInterval:numel(timepoints)-timeInterval
useDFF = dff(:, t_start:t_start+timeInterval-1);

% KS-test for the last time window
activeTag = false(size(dff, 1), 1);
for nNeuron    = 1:size(dff, 1)
    slicedDFF           = useDFF(nNeuron, :); 
    slicedDFF(isnan(slicedDFF)) = [];
    if ~isempty(slicedDFF)
        slicedDFF           = (slicedDFF - mean(slicedDFF))/std(slicedDFF);
        activeTag(nNeuron) = kstest2(-slicedDFF(slicedDFF<0), slicedDFF(slicedDFF>0), 'alpha', 0.05) && (skewness(slicedDFF)>0);
    end
end

activeNeuronMat(:, (t_start-1)/timeInterval+1) = activeTag;
end

activeTag = sum(activeNeuronMat, 2)>active_thres;


useDFF = dff(:, FA_range);
side = zeros(nCells, 1);
valid_id = ~isnan(useDFF(:, 1));
normalizedData = zscore(useDFF(valid_id, :)');
[lambda,~,~,~,~] = factoran(normalizedData,2,'maxit',1000);
colorTable = zeros(size(normalizedData, 2), 3);
colorTable(:, 1:2) = abs(lambda);
figure, hold on
for i = 1:size(normalizedData, 2)
    scatter(lambda(i, 1), lambda(i, 2), 'MarkerFaceColor', colorTable(i, :), 'MarkerEdgeColor', 'none');
end
xlabel('Factor 1');
ylabel('Factor 2');
side_valid = zeros(size(lambda, 1), 1);
side_valid(lambda(:, 1)> FA_thres) = 1;
side_valid(lambda(:, 2)> FA_thres) = 2;
side(valid_id) = side_valid;

cells_common = find(activeTag & side);
cells_active = find(activeTag & ~side);
cells_pattern = find(~activeTag & side);

save([inputFolder '.simpleTrack\active_neuron_extraction\pre_selection.mat'], 'cells_common', 'cells_active', 'cells_pattern');