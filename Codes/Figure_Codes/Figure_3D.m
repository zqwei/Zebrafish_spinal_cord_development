%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3D plot factor size as a function of cell type
%
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
addpath('../Func')
setDir;

islet_datasets = [10 13 15 16];
fs_all    = [];
islet_all = [];
mnx_all   = [];
for nFile = islet_datasets
    fileName = fileNames{nFile};
    load([tempDatDir, fileName, '.mat'], 'mnx');
    load([tempDatDir, 'Leader_', fileName, '.mat'], 'factorSize', 'islet');
    fs_all    = [fs_all; factorSize];
    islet_all = [islet_all; islet];
    mnx_all   = [mnx_all; mnx];
end

y1 = fs_all(islet_all==1 & mnx_all==1);
y2 = fs_all(islet_all==0 & mnx_all==1);
y3 = fs_all(mnx_all==0);
y2(numel(y2)+1:numel(y1)) = NaN;
y3(numel(y3)+1:numel(y1)) = NaN;
figure,
distributionPlot([y1, y2, y3], 'showMM', 0, 'globalNorm', 3, 'addSpread', 1, 'histOpt', 1, 'divFactor', 0.9);
hold on
set(gca, 'xTickLabel', {'islet+/mnx+', 'islet-/mnx+', 'mnx-'});
ylim([0 20])
set(gca, 'YTick', [1 10 20])
setPrint(8, 6, [plotDir,  'Figure_3D Leader IHC'], 'pdf');
close