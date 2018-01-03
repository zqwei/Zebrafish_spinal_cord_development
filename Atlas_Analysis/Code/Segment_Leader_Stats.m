addpath('../Func');
setDir;
fileList = [3,4,10,12,13,15,16,23];

locationList = [];
orderList = [];
for nFile = fileList
    load([TempDataDir '/tmp_' dataset{nFile} '.mat']);
    load([DirNames{nFile} '\data.mat'], 'tracks', 'timePoints');
    refTP = timePoints(find(mod(timePoints, 1200)==0, 1, 'last'));
    points = squeeze(tracks(:, refTP, :));
    actOrder = nan(size(points, 1), 1);
    metric = halfEVTime;
    for seg = floor(min(x-0.5)):ceil(max(x-0.5))
        currentSegLeft = find(y<0 & x>=seg+0.5 & x<seg+1.5 & mnx==1 & ~isnan(metric));
        [~, ~, ic] = unique(metric(currentSegLeft));
        actOrder(currentSegLeft) = ic;
        currentSegRight = find(y>0 & x>=seg+0.5 & x<seg+1.5 & mnx==1 & ~isnan(metric));
        [~, ~, ic] = unique(metric(currentSegRight));
        actOrder(currentSegRight) = ic;
    end
    locationList = [locationList; x(~isnan(actOrder))];
    orderList = [orderList; actOrder(~isnan(actOrder))];
end
locationList = locationList-0.5 - floor(locationList-0.5);
bins = 1:max(orderList)+1;
nbins = 20;
a_bins = linspace(0, 1-1/nbins, nbins);

count = polar_histogram(orderList, locationList, bins, a_bins);
figure, plot(a_bins(1:end-1), count');
legend(num2str(bins(1:end-1)'));
hold on
plot([0.5, 0.5], [0, 25], '--k');
hold off