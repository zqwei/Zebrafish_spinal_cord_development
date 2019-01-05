%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell Tracking Step 1 -- Automated cell tracking
% Update cell centroid by region growing, when merge detected, segment
% again based on k-means clustering
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%

function Cell_Tracking_s1(inputFolder, ratio)
addpath('../Func');
load([inputFolder 'signal\periods.mat'], 'nPeriod', 'periodStart');
timepoints = nPeriod-2:-1:0;
windowSize = 10;
maxDistance = 5;

%% main loop
GTFile = [inputFolder '\signal\Period_TM' num2str(timepoints(1), '%.6d') '_curated.mdf'];
inputFile = [inputFolder '\signal\Period_TM??????.klb'];
outputFolder = [inputFolder '\tracking'];

masterClock = tic;

seedPoints = readMDF(GTFile);

nTimepoints = numel(timepoints);
nCells = size(seedPoints, 1);
tracks = cell(nTimepoints, 1);

timeArray = zeros(nTimepoints, 1);
timeArray(1) = toc(masterClock);

disp(['Process time point 1/' num2str(nTimepoints) ' in ' num2str(timeArray(1), '%.2f') ' seconds']);

[tracks{1}, validLabel] = updateSeeds(recoverFilenameFromPattern(inputFile, timepoints(1)), seedPoints, ones(nCells, 1), windowSize, [], Inf);

terminationTimepoint = Inf(nCells, 1);
terminationTimepoint(validLabel==0) = 1;


for i = 2:numel(timepoints)
    timeArray(i) = toc(masterClock);
    fileName = recoverFilenameFromPattern(inputFile, timepoints(i));
    disp(['Process time point ' num2str(i) '/' num2str(numel(timepoints)) ' in ' num2str(timeArray(i) - timeArray(i-1), '%.2f') ' seconds. Valid cell number: ' num2str(sum(validLabel>0))]);
    [tracks{i}, validLabel_new] = updateSeeds(fileName, tracks{i-1}, validLabel, windowSize, ratio, maxDistance);
    terminationTimepoint(validLabel - validLabel_new > 0) = i;
    validLabel = validLabel_new;
end

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

cellNum = zeros(nTimepoints, 1);
for i = 1:nTimepoints
    cellNum(i) = sum(terminationTimepoint>=i);
end
figure, plot(cellNum);
save(fullfile(outputFolder, 'tracks.mat'), 'tracks', 'timepoints', 'GTFile', 'inputFile', 'validLabel', 'terminationTimepoint', 'timeArray');
