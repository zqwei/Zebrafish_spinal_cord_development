%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.3 Violin plot molecular identity of leaders vs. non-leaders
%
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org

function Leader_v2_3(datasets)
addpath('../Func');
setDir;

% plot distribution of cell types by neuronType
type_all = [];
fs_all = [];
pal_all = [];
for nFile = datasets;
    fileName       = fileNames{nFile};
    load([tempDatDir, 'Leader_' fileName, '.mat'], 'factorSize', 'islet', 'preActLevel');
    load([tempDatDir, fileName, '.mat'], 'mnx');
    type           = nan(size(mnx));
    type(islet==1 & mnx==1) = 1;
    type(islet==0 & mnx==1) = 2;
    type(mnx==0) = 3;
    type_all       = [type_all; type];
    fs_all         = [fs_all; factorSize;];
    pal_all        = [pal_all; preActLevel];
end

% violinPlot(fs_all, type_all, [0 20]);
% setPrint(8, 6, [plotDir,  'LeaderViolinPlot_factorSize'], 'pdf');
% close

violinPlot(pal_all, type_all, [-0.2 1.2]);
setPrint(8, 6, [plotDir,  'LeaderViolinPlot_preActLevel'], 'pdf');
% close
end

function violinPlot(me, type, range)
type(isnan(me)) = [];
me(isnan(me)) = [];

y1 = me(type==1);
y2 = me(type==2);
y3 = me(type==3);
y2(numel(y2)+1:numel(y1)) = NaN;
y3(numel(y3)+1:numel(y1)) = NaN;
figure,
distributionPlot([y1, y2, y3], 'showMM', 0, 'globalNorm', 3, 'addSpread', 1, 'histOpt', 1, 'divFactor', 0.9);
hold on
set(gca, 'xTickLabel', {'islet+', 'islet-', 'mnx-'})
ylim(range)
hold off

end