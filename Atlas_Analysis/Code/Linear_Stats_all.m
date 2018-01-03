function Linear_Stats_all()
addpath('../Func');
setDir;

%% plotting pioneer factor onto polar plot
nbins = 6;

control_datasets = [3, 4, 7, 12, 10, 11, 13, 15, 16];
MO_datasets = [17 18 19 20 21] ;

% option 1: use Ziqiang's definition of neuronType
% bins = 0.5:3;
% metric_cmd = ' me=neuronType; me(me==3)=2;';
% bins = 0.5:4;
% metric_cmd = ' me=neuronType; ';

% % option 2: use factorSize definition
bins = [0 1 1.2];
metric_cmd = ' me=exp(1-factorSize);';

% option 3: use neuronActTime
% bins = [0 0.99 1.2];
% metric_cmd = ' me=neuronActTime;';

% for control datasets
Linear_Stats_Boot_Std(control_datasets, nbins, bins, metric_cmd);
set(gcf, 'InvertHardCopy', 'off', 'PaperPositionMode', 'Auto');   
% print([PlotDir '/seg_distr_dkbg/WT_segRatio_leader.pdf'], '-dpdf', '-r0');

% for MO
Linear_Stats_Boot_Std(MO_datasets, nbins, bins, metric_cmd);
% print([PlotDir '/seg_distr_dkbg/MO_segRatio_leader.pdf'], '-dpdf', '-r0');


%% plotting individual circular statistics

% file_set = [1 2 3 4 5 8 9];
% r_mean = zeros(7, 2);
% r_dist = zeros(7, 2);
% for i = 1:numel(file_set);
%     nFile = file_set(i);
%     load([TempDataDir '/tmp_' dataset{nFile} '.mat']);
%     me = exp(1-factorSize);
%
%     me_plot = [me(x>=1 & x<=size(ra, 1) & ~isnan(me))];
%     x_plot = [x(x>=1 & x<=size(ra, 1) & ~isnan(me))];
%     rad = (x_plot-floor(x_plot))*2*pi;
%     r_mean(i, 1) = circ_mean(rad(me_plot<thres));
%     r_mean(i, 2) = circ_mean(rad(me_plot>=thres));
%     r_dist(i, 1) = circ_r(rad(me_plot<thres));
%     r_dist(i, 2) = circ_r(rad(me_plot>=thres));
% end



% plot distribution of cell types by neuronType
islet_datasets = [10 13 15 16];
islet_all = [];
mnx_all = [];
me_all = [];
for nFile = islet_datasets;
    load([TempDataDir '/tmp_' dataset{nFile} '.mat'], 'factorSize', 'islet', 'mnx', 'x');
    load([DirNames{nFile} '/LONOLoading_v_0_1.mat'], 'neuronType');
    load(['../../NetworkDynamicsTempDat/NeuronActTime' dataset{nFile} '.mat'], 'neuronActTime');
    eval(metric_cmd);
    
    if nFile == 16
        me(x<1) = NaN;
    end
    islet_all = [islet_all; islet];
    mnx_all = [mnx_all; mnx];
    me_all = [me_all; me;];
    
%     if sum(me>0.7 & mnx==0)
%         disp([dataset{nFile} ' cell ' num2str(find(me>0.7 & mnx==0))]);
%     end
    
    
end

islet_all(isnan(me_all)) = [];
mnx_all(isnan(me_all)) = [];
me_all(isnan(me_all)) = [];

count = zeros(numel(bins)-1, 3);
for type = 1: numel(bins) - 1
    %     select = me_all >= bins(type) & me_all<bins(type+1);
    select = me_all >= bins(type);
    count(type, 1) = sum(select & mnx_all==1 & islet_all==1);
    count(type, 2) = sum(select & mnx_all==1 & islet_all==0);
    count(type, 3) = sum(select & mnx_all==0);
end
figure,
set(0, 'defaultaxeslayer', 'top')
whitebg('whilte');
set(gcf, 'Color', 'black');
% set(0 ,'defaultfigurecolor', 'black')
bar(count,'stacked');
set(gca, 'XTickLabel', {'type 1', 'type 2', 'type3'});
legend({'islet+', 'islet-', 'mnx-'});
title('activation type and cell type');
ylabel('count');
% export_fig([PlotDir '/seg_distr_dkbg/IHC_leader.pdf'], '-nocrop');
% Version for normalization
% figure,
% bar(count./repmat(sum(count, 2), 1, 3), 'stacked');
% set(gca, 'XTickLabel', {'type 1', 'type 2', 'type3'});
% legend({'islet+', 'islet-', 'mnx-'});
% title('activation type and cell type');
% ylabel('normalized percentage');

me_all = 1-log(me_all);
y1 = me_all(islet_all==1 & mnx_all==1);
y2 = me_all(islet_all==0 & mnx_all==1);
y1(numel(y1)+1:numel(y2)) = NaN;
y2(numel(y2)+1:numel(y1)) = NaN;
figure,
distributionPlot([y1, y2], 'showMM', 0, 'globalNorm', 3, 'addSpread', 1, 'histOpt', 1, 'divFactor', 0.9);
hold on
set(gca, 'xTickLabel', {'islet+', 'islet-'})
ylim([0 20])
set(gca, 'YTick', [1 10 20])
% export_fig([PlotDir '/islet_fs_dist_' dataset{nFile} '.tif'], '-nocrop');
hold off
% close

end

function Linear_Stats_Boot_Std(datasets, nbins, bins, metric_cmd)
addpath('../Func');
setDir;
a_bins = linspace(0, 1, nbins+1);
leg = linspace(0, 360 - 360/nbins, nbins);
count = zeros(numel(datasets), nbins, numel(bins)-1); %nFile x nSector x nType
for i = 1:numel(datasets)
    nFile = datasets(i);
    load([TempDataDir '/tmp_' dataset{nFile} '.mat'], 'x', 'factorSize');
    load([DirNames{nFile} '/LONOLoading_v_0_1.mat'], 'neuronType');
    load(['../../NetworkDynamicsTempDat/NeuronActTime' dataset{nFile} '.mat'], 'neuronActTime');
    eval(metric_cmd);
    flag_plot = x>=1 & x<=floor(max(x)) & ~isnan(me);
    me = me(flag_plot);
    x = x(flag_plot);
    count(i, :, :) = polar_histogram(me, x-floor(x), bins, a_bins);
end

% bootstraping std
bootStd  = zeros(nbins, numel(bins)-1);
bootMean = zeros(nbins, numel(bins)-1);
for nType = 1:numel(bins)-1
    bootSum = bootstrp(1000, @sum, squeeze(count(:, :, nType)));
    bootStd(:, nType) = std(bootSum, [], 1);
    bootMean(:, nType) = mean(bootSum);
end

count_all = squeeze(sum(count, 1));
bootStd = bootStd./repmat(sum(count_all, 1), nbins, 1);
count_all = count_all./repmat(sum(count_all, 1), nbins, 1);
%     count_plot = [count_all-bootStd, count_all, count_all+bootStd];
%     spider(count_plot, '', repmat([0, max(count_plot(:))], nbins, 1), strtrim(cellstr(num2str(leg'))'));
figure, 
set(0, 'defaultaxeslayer', 'top')
whitebg('white');
set(gcf, 'Color', 'black');
% set(0 ,'defaultfigurecolor', 'black')
hold on
anglebins = linspace(0, 360, nbins+1);
index = circshift((1:nbins)', nbins/2);
index = [index; index(1)];
cmap = lines(2);
hl = zeros(size(count_all, 2), 1);
for i = 1:size(count_all, 2)
    h = area(anglebins, [count_all(index, i)-bootStd(index, i), 2*bootStd(index, i)], 'EdgeColor', 'none');
    set(h(2), 'FaceColor', cmap(i, :), 'FaceAlpha', 0.5);
    set(h(1), 'FaceColor', 'k', 'FaceAlpha', 0.5);
    hl(i) = plot(anglebins, count_all(index, i), 'Color', cmap(i, :));
end
set(gca, 'XTick', anglebins, 'XTickLabel', strtrim(cellstr(num2str(linspace(-180, 180, nbins+1)'))));
xlim([0, 360]);
legend(hl, {'non-leader', 'leader'})
%     % mean and regular std - not working because of negative values!!
%     count_all = squeeze(mean(count, 1));
%     count_plot = [count_all-squeeze(std(count)), count_all, count_all+squeeze(std(count))];
%     spider(count_plot, '', repmat([0, max(count_plot(:))], nbins, 1), strtrim(cellstr(num2str(leg'))'));
end