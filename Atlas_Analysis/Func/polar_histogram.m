function polar_histogram(stats, t, c_bins, t_bins)
% Function to make polar histogram of stats(t), where stats is categorized
% to stats_bins in different color, and shown in different sectors
% according to the distribution of t in nbins
% By Yinan Wan
% input paras:
% stats: input metrics of the input stats, N x 1 vector
% t: input metrics of the polar angle, N x 1 vector
% c_bins: bins to make histograms in different colors
% t_bins: bins to make histograms in different ploar sectors


leg = linspace(0, 360 - 360/(numel(t_bins) - 1), numel(t_bins)-1);

count = zeros(numel(t_bins)-1, numel(c_bins)-1);

for i = 1:numel(t_bins)-1
    select = find(t>=t_bins(i) & t<t_bins(i+1));
    if ~isempty(select)
        count(i, :) = histcounts(stats(select), c_bins);
    end
end

% if normalization is used
% count = count./repmat(max(count), numel(t_bins)-1, 1);

spider(count, '', repmat([0, max(count(:))], numel(t_bins)-1, 1), strtrim(cellstr(num2str(leg'))'));

end