addpath('../Func');
setDir;

%% plotting pioneer factor onto polar plot

me_plot = [];
x_plot = [];
for nFile = [3, 4, 7, 12, 10, 11, 13, 15, 16]
% for nFile =  [17 18 19 20 21] 

    load([TempDataDir '/tmp_' dataset{nFile} '.mat']);
    me = exp(1-factorSize);

    me_plot = [me_plot; me(x>=1 & x<=floor(max(x)) & ~isnan(me))];
    x_plot = [x_plot; x(x>=1 & x<=floor(max(x)) & ~isnan(me))];
end


minR = 0;
maxR = 1;
ticks = 0.2:0.2:1;
thres = exp(-1);
nbins = 12;

% spider plot
bins = [0 thres 1];
a_bins = linspace(0, 1, nbins+1);
count = zeros(numel(a_bins)-1, numel(bins)-1);
for i = 1:numel(a_bins)-1
    select = find(x_plot-floor(x_plot)>=a_bins(i) & x_plot-floor(x_plot)<a_bins(i+1));
    if ~isempty(select)
        count(i, :) = histcounts(me_plot(select), bins);
    end
end
leg = 0:30:330;
spider(count, '', repmat([0, max(count(:))], numel(a_bins)-1, 1), strtrim(cellstr(num2str(leg'))'));
% spider(count(:, 1), '', repmat([0, max(count(:,1))], numel(a_bins)-1, 1), strtrim(cellstr(num2str(leg'))'));
% spider(count(:, 2), '', repmat([0, max(count(:,2))], numel(a_bins)-1, 1), strtrim(cellstr(num2str(leg'))'));

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



% %% plotting islet with factorSize
% for nFile = [10 13 15 16]
% islet_all = [];
% mnx_all = [];
% me = [];
% figure,
% clf('reset')
% % for nFile = [10 13 15 16]
%     load([TempDataDir '/tmp_' dataset{nFile} '.mat']);
%     islet_all = [islet_all; islet];
%     mnx_all = [mnx_all; mnx];
%     me = [me; factorSize;];
% % end
% y1 = me(islet_all==1 & mnx_all==1);
% y2 = me(islet_all==0 & mnx_all==1);
% y1(numel(y1)+1:numel(y2)) = NaN;
% y2(numel(y2)+1:numel(y1)) = NaN;
% distributionPlot([y1, y2], 'showMM', 0, 'globalNorm', 3, 'addSpread', 1);
% hold on
% % scatter(ones(numel(y1), 1)+randn(numel(y1), 1)/10, y1,10, 'ok', 'MarkerFaceColor', [0.5, 0.5, 0.5]);
% % scatter(ones(numel(y2), 1)+randn(numel(y2), 1)/10+1, y2, 10, 'ok', 'MarkerFaceColor', [0.5, 0.5, 0.5]);
% set(gca, 'xTickLabel', {'islet+', 'islet-'})
% % ylim([0 10])
% % export_fig([PlotDir '/islet_dist_' dataset{nFile} '.tif'], '-nocrop');
% % hold off
% % close
% end