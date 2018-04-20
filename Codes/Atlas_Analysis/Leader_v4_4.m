%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.4 
% Monte-Carlo simulation
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
function Leader_v4_4(control_datasets)

addpath('../Func')
setDir;

preActAll = [];
for i = 1:numel(control_datasets)
    nFile = control_datasets(i);
    fileName = fileNames{nFile};
    load([tempDatDir, 'Leader_', fileName, '.mat'], 'preActLevel');
            preActAll = [preActAll; preActLevel];
end

preActAll(isnan(preActAll)) = [];
bins = 0:0.01:1;
[f, xi] = ksdensity(preActAll, 'bandwidth', 0.02);
percentage = histcounts(preActAll, bins)/numel(preActAll);
figure, bar(bins(1:end-1), percentage/max(percentage), 'histc');
[~, lct] = findpeaks(-f);
disp(num2str(xi(lct)))
hold on,
plot(xi, f/max(f), 'r');
plot([xi(lct), xi(lct)], [0, 1], 'k--')
xlabel('activity level before factored');
ylabel('percentage');
xlim([0, 1])
setPrint(8, 6, [plotDir,  'PreActLevelDistr_all'], 'pdf');
% close;
