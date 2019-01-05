%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Signal Extraction Step 5: Extract signal from curated tracks
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
function Signal_Extraction_s5(inputFolder, ratio3D)
addpath('../Func');
dataFile = [inputFolder '.simpleTrack\tracking\tracks_smoothed.mat'];
periodFile = [inputFolder '.simpleTrack\signal\periods.mat'];
load(periodFile); 
load(dataFile); % load 4D matrix tracks_smoothed of nLineage x nPeriod x 3
inputFile = [inputFolder '\SPM00\TM??????\SPM00_TM??????_CM00_CHN00.klb'];
outputFolder = [inputFolder '.simpleTrack\profile'];
localRun = [1 12]; % localRun(1) 0: use the cluster 1: run on local machine (specify number of pool workers in localRun(2))

radius = 4;
jobMemory = 1;


% interpolate to get the full range of cell tracks
timepoints_full = min(periodStart(timePoints+1)):max(periodStart(timePoints+1));
nTimepoints = numel(timepoints_full);
nCells = size(tracks_smoothed, 1);
tracks_all = nan(nCells, nTimepoints, 3);

for i = 1:nCells
    tracks_all(i, :, 1) = interp1(periodStart(timePoints+1), tracks_smoothed(i, :, 1), timepoints_full);
    tracks_all(i, :, 2) = interp1(periodStart(timePoints+1), tracks_smoothed(i, :, 2), timepoints_full);
    tracks_all(i, :, 3) = interp1(periodStart(timePoints+1), tracks_smoothed(i, :, 3), timepoints_full);
end

% calculate sphere coordinates based on meshgrid (on the isotropic scale)
stackDimensions = size(readImage(recoverFilenameFromPattern(inputFile, timepoints_full(1))));
[rr, cc] = meshgrid(-radius:radius, -radius:radius);
[cmx, cmy] = ind2sub([2*radius+1, 2*radius+1], find(rr.^2 + cc.^2 <= radius^2));



if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Ignore existing time point
tagExist = false(nTimepoints, 1);
for i = 1:nTimepoints
    if exist([outputFolder '\profile.TM' num2str(timepoints_full(i), '%.6d') '.mat'], 'file')
        tagExist(i) = 1;
    end
end
timepoints_full = timepoints_full(tagExist==0);
nTimepoints = numel(timepoints_full);
tracks_all = tracks_all(:, tagExist==0, :);

% rename "tracks_all" and "timepoints_full" to make sure we can use extractSignal_fast
tracks_smoothed = tracks_all;
timepoints = timepoints_full;

currentTime = clock;
timeString = [...
    num2str(currentTime(1)) num2str(currentTime(2), '%.2d') num2str(currentTime(3), '%.2d') ...
    '_' num2str(currentTime(4), '%.2d') num2str(currentTime(5), '%.2d') num2str(round(currentTime(6) * 1000), '%.5d')];
parameterDatabase = [pwd '\jobParameters.extractSignal_all' timeString '.mat'];

save(parameterDatabase, 'timepoints', 'inputFile', 'radius', 'tracks_smoothed', 'nCells', 'cmx', 'cmy', 'outputFolder', 'nTimepoints', 'stackDimensions', 'ratio3D');

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
collectResult(parameterDatabase)



