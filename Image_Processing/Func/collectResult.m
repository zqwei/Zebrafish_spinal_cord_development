%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% collectResult.m
% Collect result from parallel processing
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%

function collectResult(databasename)
load(databasename);
profile_all = zeros(nCells, nTimepoints);
for i = 1:nTimepoints
    load([outputFolder '\profile.TM' num2str(timepoints(i), '%.6d') '.mat']);
    disp(['loading time point TM' num2str(timepoints(i), '%.6d')]);
    profile_all(:, i) = profile;
end
save(fullfile(outputFolder, 'profile.mat'), 'profile_all', 'timepoints', 'tracks_smoothed', 'nCells');
for i = 1:nTimepoints
    delete([outputFolder '\profile.TM' num2str(timepoints(i), '%.6d') '.mat']);
end
