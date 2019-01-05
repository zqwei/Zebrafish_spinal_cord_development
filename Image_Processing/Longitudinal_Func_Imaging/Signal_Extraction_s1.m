%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Signal Extraction Step 1: Extraction of automated tracking result for pre-selection
% Use a time_range of time points from the end of recording
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%
function Signal_Extraction_s1(inputFolder, time_range, ratio3D)
addpath('../Func');
dataFile = [inputFolder '.simpleTrack\tracking\tracks.mat'];
periodFile = [inputFolder '.simpleTrack\signal\periods.mat'];
load(periodFile);
load(dataFile);

inputFile = [inputFolder '\SPM00\TM??????\SPM00_TM??????_CM00_CHN00.klb'];
outputFolder = [inputFolder '.simpleTrack\active_neuron_extraction'];
localRun = [1 12]; % localRun(1) 0: use the cluster 1: run on local machine (specify number of pool workers in localRun(2))

radius = 4;
jobMemory = 1;

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% use all time points from period 1 to preiod nPeriod-1
timepoints = periodStart(1): periodStart(end-1)- 1;
nTimepoints = numel(timepoints);

nCells = size(tracks{1}, 1);
tracks_smoothed = zeros(nCells, nTimepoints, 3);
for i = 1:nCells
    currentCellTrack = zeros(nPeriod-1, 3);
    for j = 1:nPeriod-1
        if j < terminationTimepoint(i)
            currentCellTrack(nPeriod-j, :) = tracks{j}(i, :);
        else
            currentCellTrack(nPeriod-j, :) = NaN;
        end
    end
    tracks_smoothed(i, :, :) = interp1(periodStart(1:end-1), currentCellTrack, timepoints);
end


stackDimensions = size(readImage(recoverFilenameFromPattern(inputFile, timepoints(1))));
[rr, cc] = meshgrid(-radius:radius, -radius:radius);
[cmx, cmy] = ind2sub([2*radius+1, 2*radius+1], find(rr.^2 + cc.^2 <= radius^2));

% Ignore existing time point - IMPORTANT: use here to only calculate the
% last 7200 time points
tagExist = false(nTimepoints, 1);
tagExist(1:end-time_range) = 1;
for i = 1:nTimepoints
    if exist([outputFolder '\profile.TM' num2str(timepoints(i), '%.6d') '.mat'], 'file')
        tagExist(i) = 1;
    end
end
timepoints = timepoints(tagExist==0);
nTimepoints = numel(timepoints);
tracks_smoothed = tracks_smoothed(:, tagExist==0, :);

currentTime = clock;
timeString = [...
    num2str(currentTime(1)) num2str(currentTime(2), '%.2d') num2str(currentTime(3), '%.2d') ...
    '_' num2str(currentTime(4), '%.2d') num2str(currentTime(5), '%.2d') num2str(round(currentTime(6) * 1000), '%.5d')];
parameterDatabase = [pwd '\jobParameters.extractSignal_preselect' timeString '.mat'];

save(parameterDatabase, 'timepoints', 'inputFile', 'radius', 'tracks_smoothed', 'nCells', 'cmx', 'cmy','outputFolder', 'nTimepoints', 'stackDimensions', 'ratio3D');
tic,
if localRun(1)
    if matlabpool('size') > 0
        matlabpool('close');
    end;
    matlabpool(localRun(2));
    
    parfor i = 1:nTimepoints
        disp(['Processing time point TM' num2str(timepoints(i), '%.6d') ' on local machine']);
        extractSignal_fast(parameterDatabase, i);
    end
    matlabpool('close');
else
    disp(['Submitting parametric cluster job for ' num2str(nTimepoints) ' time point(s).']);
    cmdFunction = ['extractSignal_fast(''' parameterDatabase ''', *)'];
    cmd = ['job submit /scheduler:keller-cluster.janelia.priv /user:simview ' ...
        '/parametric:1-' num2str(nTimepoints) ':1 runMatlabJob.cmd """' pwd '""" """' cmdFunction '"""'];
    [status, systemOutput] = system(cmd);
    disp(['System response: ' systemOutput]);
end
toc
collectResult(parameterDatabase)


