%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.2 Plot molecular identity of leaders vs. non-leaders
%
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org

function Leader_v2_2(datasets, preActThres)
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
    fs_all         = [fs_all; exp(1-factorSize);];
    pal_all        = [pal_all; preActLevel];
end

% leaderBarPlot(fs_all, exp(1-fsThres), type_all, 0);
% setPrint(8*2, 6, [plotDir,  'LeaderBarPlot_factorSize_thres_' num2str(fsThres)], 'pdf');
% % close;
leaderBarPlot(pal_all, preActThres, type_all, 1);
setPrint(8*2, 6, [plotDir,  'LeaderBarPlot_preActLevel_thres_' num2str(preActThres)], 'pdf');
% close;
end
function leaderBarPlot(me, meThres, type, normTag)
type(isnan(me)) = [];
me(isnan(me)) = [];
leader = me>meThres; % 1:leader, 0:non-leader

% (leader, non-leader) with stacked type
count = zeros(2, 3);
for nType = 1:3
    count(1, nType) = sum(leader & type==nType);
    count(2, nType) = sum(~leader & type==nType);
end
figure,
subplot(1, 2, 1)
count_plot = count;
if normTag
    count_plot = count./repmat(sum(count, 2), 1, 3);
end
bar(count_plot,'stacked');
set(gca, 'XTickLabel', {'leader', 'non-leader'});
legend({'islet+', 'islet-', 'mnx-'});
ylabel('count');

subplot(1, 2, 2)
count_plot = count;
if normTag
    count_plot = count./repmat(sum(count), 2, 1);
end
bar(count_plot','stacked');
set(gca, 'XTickLabel', {'islet+', 'islet-', 'mnx-'});
legend({'leader', 'non-leader'});
ylabel('count');
end