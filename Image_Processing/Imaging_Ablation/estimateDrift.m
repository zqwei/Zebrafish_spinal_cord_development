function estimateDrift(nFile)
addpath('Func');
setDir;
timeStep = 120;
smoothWindow = 60;
refTP = [0:timeStep:inputFolderList{nFile, 2}, inputFolderList{nFile, 2}];
drift = zeros(numel(refTP), 3);
for t = 1:numel(refTP)-1
    M = readTrsf([inputFolderList{nFile, 1} '/DriftCorrection/affine_' num2str(refTP(t), '%.6d') '.trsf']);
    drift(t, :) = M(1:3, 4)';
end
drift(:, 3) = drift(:, 3)/6;
fullTP = 0: inputFolderList{nFile, 2};
driftSmoothed = zeros(numel(fullTP), 3);
for dim = 1:3
    driftSmoothed(:, dim) = interp1(refTP, drift(:, dim), fullTP);
    driftSmoothed(:, dim) = smooth(driftSmoothed(:, dim), smoothWindow);
end
figure, hold on
plot(refTP, drift, 'ok');
h = plot(fullTP, driftSmoothed);
legend(h, {'x', 'y', 'z'});
hold off
title(['fish ' num2str(nFile) ' before']);
fullDrift.before.TP = fullTP;
fullDrift.before.drift = driftSmoothed;

refTP = 0:timeStep:inputFolderList{nFile, 4};
drift = zeros(numel(refTP), 3);
for t = 1:numel(refTP)
    M = readTrsf([inputFolderList{nFile, 3} '/DriftCorrection/affine_' num2str(refTP(t), '%.6d') '.trsf']);
    drift(t, :) = M(1:3, 4)';
end
drift(:, 3) = drift(:, 3)/6;
fullTP = 0:inputFolderList{nFile, 4};
driftSmoothed = zeros(numel(fullTP), 3);
for dim = 1:3
    driftSmoothed(:, dim) = interp1(refTP, drift(:, dim), fullTP);
    driftSmoothed(:, dim) = smooth(driftSmoothed(:, dim), smoothWindow);
end
figure, hold on
plot(refTP, drift, 'ok');
h = plot(fullTP, driftSmoothed);
hold off
legend(h, {'x', 'y', 'z'});
title(['fish ' num2str(nFile) ' after']);
fullDrift.after.TP = fullTP;
fullDrift.after.drift = driftSmoothed;
save([inputFolderList{nFile, 1} '/DriftCorrection/estimatedDrift.mat'], 'fullDrift');