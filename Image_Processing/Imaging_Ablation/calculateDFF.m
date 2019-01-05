function calculateDFF(nFile)
addpath('Func');
setDir;
baselineWindowSize = 40;
baselinePrc = 20;
background = 90;

folderIDList = [1, 3];

for folderID = 1:numel(folderIDList)
    nFolder = folderIDList(folderID);
    load([inputFolderList{nFile, nFolder} '/Profile/profile.mat']);
    baseline = zeros(size(profile_all));
    for i = 1:nCells
        baseline(i, :) = running_percentile(profile_all(i, :), baselineWindowSize, baselinePrc);
    end
    dff = (profile_all-baseline)./(baseline - background);
    save([inputFolderList{nFile, nFolder} '/Profile/profile.mat'], 'dff', '-append');
end