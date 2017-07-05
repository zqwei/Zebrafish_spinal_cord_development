addpath('../Func');
setDir;

%% plotting pioneer factor onto polar plot 
nbins = 8;
a_bins = linspace(0, 1, nbins+1);
control_datasets = [3, 4, 7, 12, 10, 11, 13, 15, 16];
MO_datasets = [17 18 19 20 21] ;

% option 1: use Ziqiang's definition of neuronType
% bins = 0.5:3;
% metric_cmd = ' me=neuronType; me(me==2)=1; me(me==3)=2;';
bins = 0.5:4;
metric_cmd = ' me=neuronType; ';


% option 2: use factorSize definition
% bins = [0 exp(-1) 1];
% metric_cmd = ' me=exp(1-factorSize);';


% for control datasets
me_plot = [];
x_plot = [];
for nFile = control_datasets
    load([TempDataDir '/tmp_' dataset{nFile} '.mat'], 'x', 'factorSize');
    load([DirNames{nFile} '/LONOLoading_v_0_1.mat'], 'neuronType');
    eval(metric_cmd);
    me_plot = [me_plot; me(x>=1 & x<=floor(max(x)) & ~isnan(me))];
    x_plot = [x_plot; x(x>=1 & x<=floor(max(x)) & ~isnan(me))];
end
polar_histogram(me_plot, x_plot-floor(x_plot), bins, a_bins);

% for MO
me_plot = [];
x_plot = [];
for nFile =  MO_datasets
    load([TempDataDir '/tmp_' dataset{nFile} '.mat'], 'x', 'factorSize');
    load([DirNames{nFile} '/LONOLoading_v_0_1.mat'], 'neuronType');
    eval(metric_cmd);
    me_plot = [me_plot; me(x>=1 & x<=floor(max(x)) & ~isnan(me))];
    x_plot = [x_plot; x(x>=1 & x<=floor(max(x)) & ~isnan(me))];
end

polar_histogram(me_plot, x_plot-floor(x_plot), bins, a_bins);

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
    load([TempDataDir '/tmp_' dataset{nFile} '.mat'], 'factorSize', 'islet', 'mnx');
    load([DirNames{nFile} '/LONOLoading_v_0_1.mat'], 'neuronType');
    eval(metric_cmd);
    islet_all = [islet_all; islet];
    mnx_all = [mnx_all; mnx];
    me_all = [me_all; me;];
end

islet_all(isnan(me_all)) = [];
mnx_all(isnan(me_all)) = [];
me_all(isnan(me_all)) = [];

count = zeros(numel(bins)-1, 3);
for type = 1: numel(bins) - 1
    select = me_all >= bins(type) & me_all<bins(type+1);
    count(type, 1) = sum(me_all(select & mnx_all==1 & islet_all==1));
    count(type, 2) = sum(me_all(select & mnx_all==1 & islet_all==0));
    count(type, 3) = sum(me_all(select & mnx_all==0));
end
figure, 
bar(count,'stacked');
set(gca, 'XTickLabel', {'type 1', 'type 2', 'type3'});
legend({'islet+', 'islet-', 'mnx-'});
title('activation type and cell type');
ylabel('count');

figure, 
bar(count./repmat(sum(count, 2), 1, 3), 'stacked');
set(gca, 'XTickLabel', {'type 1', 'type 2', 'type3'});
legend({'islet+', 'islet-', 'mnx-'});
title('activation type and cell type');
ylabel('normalized percentage');