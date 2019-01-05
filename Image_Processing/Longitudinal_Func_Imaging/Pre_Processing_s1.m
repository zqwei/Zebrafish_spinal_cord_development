%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre_Processing_v1 -- Downsampling of functional imaging data
% Max-intensity projection over time
% to run with parallel processing in cluster environment
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%

function Pre_Processing_s1(inputFilePattern, outputFolder, timepoints, timeWindow)
addpath('../Func');

localRun = [1 12];

%% main loop
outputFolder = [outputFolder '\signal'];
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

disp(' ');
disp('Setting up processing environment');

masterClock = tic;

periodStart = 0 : timeWindow : timepoints(end);
nPeriod = numel(periodStart);

currentTime = clock;
timeString = [...
    num2str(currentTime(1)) num2str(currentTime(2), '%.2d') num2str(currentTime(3), '%.2d') ...
    '_' num2str(currentTime(4), '%.2d') num2str(currentTime(5), '%.2d') num2str(round(currentTime(6) * 1000), '%.5d')];
parameterDatabase = [pwd '\jobParameters.projectStack.' timeString '.mat'];

save(parameterDatabase, 'timepoints', 'inputFilePattern', 'periodStart', 'nPeriod', 'outputFolder', 'timeWindow');
tic,

validQueue = ones(nPeriod, 1);
for i = 1:nPeriod
    if exist([outputFolder '\Period_TM' num2str(i-1, '%.6d') ',klb'], 'file');
            validQueue(i) = 0;
    end
end

jobQueue = find(validQueue);
if localRun(1)
    if matlabpool('size') > 0
        matlabpool('close');
    end;
    matlabpool(localRun(2));
    parfor i = 1:numel(jobQueue)
        projectStack(parameterDatabase, jobQueue(i));
    end
else
    disp(['Submitting parametric cluster job for ' num2str(nPeriod) ' time point(s).']);
    cmdFunction = ['projectStack(''' parameterDatabase ''', *)'];
    cmd = ['job submit /scheduler:keller-cluster.janelia.priv /user:simview ' ...
        '/parametric:1-' num2str(nPeriod) ':1 runMatlabJob.cmd """' pwd '""" """' cmdFunction '"""'];
    [status, systemOutput] = system(cmd);
    disp(['System response: ' systemOutput]);
end

save([outputFolder '\periods.mat'], 'nPeriod', 'periodStart');
