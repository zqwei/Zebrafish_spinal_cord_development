%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.1 Plot linear distribution of leaders locations - based on preActLevel
%
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org

function Leader_v2_1(control_datasets, preActThres)
addpath('../Func');
setDir;

%% plotting pioneer factor onto polar plot
nbins = 6;
bins = [0 preActThres 1];

% for control datasets
Linear_Stats_Boot_Std(control_datasets, nbins, bins);
set(gcf, 'InvertHardCopy', 'off', 'PaperPositionMode', 'Auto');   
print([plotDir 'WT_segRatio_leader_preActLevel.pdf'], '-dpdf', '-r0');



end

function Linear_Stats_Boot_Std(datasets, nbins, bins)
addpath('../Func');
setDir;
a_bins = linspace(0, 1, nbins+1);
leg = linspace(0, 360 - 360/nbins, nbins);
count = zeros(numel(datasets), nbins, numel(bins)-1); %nFile x nSector x nType
x_all = [];
me_all = [];
fishNumber = [];
for i = 1:numel(datasets)
    fileName       = fileNames{datasets(i)};
    load([tempDatDir, 'Leader_' fileName, '.mat'], 'preActLevel');
    load([tempDatDir, fileName, '.mat'], 'new_x');
    x              = new_x;
    me             = preActLevel;

    flag_plot      = x>=1 & x<=floor(max(x)) & ~isnan(me);
    me             = me(flag_plot);
    x              = x(flag_plot);
    count(i, :, :) = polar_histogram(me, x-floor(x), bins, a_bins);
    x_all = [x_all; x-floor(x)];
    me_all = [me_all; me];
    fishNumber = [fishNumber; ones(numel(me), 1)*i];
end

% bootstraping std
bootStd  = zeros(nbins, numel(bins)-1);
bootMean = zeros(nbins, numel(bins)-1);
for nType = 1:numel(bins)-1
    bootSum = bootstrp(1000, @sum, squeeze(count(:, :, nType)));
    bootStd(:, nType) = std(bootSum, [], 1);
    bootMean(:, nType) = mean(bootSum);
end
pd = makedist('uniform');

[h, p] = kstest2(x_all(me_all<=bins(2)), x_all(me_all>bins(2)));
[hN, pN] = kstest(x_all(me_all<=bins(2)), 'cdf', pd);
[hL, pL] = kstest(x_all(me_all>bins(2)), 'cdf', pd);

disp(['p=', num2str(p) ', pL=', num2str(pL) ', pN=', num2str(pN)]);

hFish = zeros(numel(datasets), 3); % L vs. N, N vs. uniform, L vs. uniform
pFish = zeros(numel(datasets), 3);
for i = 1:numel(datasets)
    [hFish(i, 1), pFish(i, 1)] = kstest2(x_all(me_all <= bins(2) & fishNumber==i), x_all(me_all>bins(2) & fishNumber==i));
    [hFish(i, 2), pFish(i, 2)] = kstest( x_all(me_all <= bins(2) & fishNumber==i), 'cdf', pd);
    [hFish(i, 3), pFish(i, 3)] = kstest( x_all(me_all >  bins(2) & fishNumber==i), 'cdf', pd);
end
pFishSorted = [sort(pFish(:, 1), 'descend'), sort(pFish(:, 2), 'descend'), sort(pFish(:, 3), 'descend')];
pMax = zeros(1, 3);
for i =1:numel(datasets)
    pMax = max(pMax, pFishSorted(i, :)*i);
end
disp(['Corrected p=' num2str(pMax)]);
count_all = squeeze(sum(count, 1));
bootStd = bootStd./repmat(sum(count_all, 1), nbins, 1);
count_all = count_all./repmat(sum(count_all, 1), nbins, 1);
%     count_plot = [count_all-bootStd, count_all, count_all+bootStd];
%     spider(count_plot, '', repmat([0, max(count_plot(:))], nbins, 1), strtrim(cellstr(num2str(leg'))'));
figure, 
set(0, 'defaultaxeslayer', 'top')
whitebg('white');
% set(gcf, 'Color', 'black');
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
    set(h(1), 'FaceColor', 'w', 'FaceAlpha', 0.5);
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