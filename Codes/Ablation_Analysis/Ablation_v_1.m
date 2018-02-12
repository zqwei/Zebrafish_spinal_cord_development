%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.  evaluate the acitavion consistency (after), plot with regard to
%
% 1) anatomical location
% 2) whether it is an active cell before
% 3) decide the threshold for activation level based on the histogram above
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%

function Ablation_v_1(fishIndex)
addpath('../Func');
setDir;

fileListOffset = 24;
for nExp = 1:2
    nFile = fileListOffset + (fishIndex-1) *2 + nExp;
    fileDirName   = fileDirNames{nFile}; %#ok<*USENS>
    fileName      = fileNames{nFile};
    
    dirImageData  = [fileDirName '/'];
    load([dirImageData, 'profile.mat'], 'segAblation', 'x', 'y', 'z');
    
    if nExp == 1
        load([tempDatDir, fileName, '.mat'], 'activeNeuronMat');
        activeTagBefore =  activeNeuronMat;
    else
        copyfile([tempDatDir, fileName, '.mat'], [tempDatDir, fileName, '_tmp', '.mat']);
        activeNeuronMat = Neuron_selection_v0_short_win(nFile);
        activeLevelAfter = sum(activeNeuronMat, 2)/size(activeNeuronMat, 2);
%         Neuron_selection_v2(nFile);
    end
end


y = y/2;
z = z/max(z) * 1.8;

figure,
subplot(3, 1, 1);
hold on
scatter(x(~activeTagBefore), y(~activeTagBefore), [], activeLevelAfter(~activeTagBefore), 'filled', 'd');
scatter(x(activeTagBefore), y(activeTagBefore), 40, activeLevelAfter(activeTagBefore), 'filled', 'o', 'MarkerEdgeColor', 'r');

for i = 1:ceil(max(x))
    if ismember(i, segAblation)
        plot([i, i], [-2, 2], '-k');
    else
        plot([i, i], [-2, 2], '--k');
    end
end
xlim([1 ceil(max(x))]);
ylim([-1 1]);
colorbar
hold off

subplot(3, 1, 2);
hold on
scatter(x(~activeTagBefore & y<0), z(~activeTagBefore & y<0), [], activeLevelAfter(~activeTagBefore & y<0), 'filled', 'd');
scatter(x(activeTagBefore & y<0), z(activeTagBefore & y<0), 40, activeLevelAfter(activeTagBefore & y<0), 'filled', 'o', 'MarkerEdgeColor', 'r');

for i = 1:ceil(max(x))
    if ismember(i, segAblation)
        plot([i, i], [-2, 2], '-k');
    else
        plot([i, i], [-2, 2], '--k');
    end
end
xlim([1 ceil(max(x))]);
ylim([-0.5, 2]);
colorbar
hold off

subplot(3, 1, 3);
hold on
scatter(x(~activeTagBefore & y>0), z(~activeTagBefore & y>0), [], activeLevelAfter(~activeTagBefore & y>0), 'filled', 'd');
scatter(x(activeTagBefore & y>0), z(activeTagBefore & y>0), 40, activeLevelAfter(activeTagBefore & y>0), 'filled', 'o', 'MarkerEdgeColor', 'r');

for i = 1:ceil(max(x))
    if ismember(i, segAblation)
        plot([i, i], [-2, 2], '-k');
    else
        plot([i, i], [-2, 2], '--k');
    end
end
xlim([1 ceil(max(x))]);
ylim([-0.5, 2]);
colorbar
hold off
export_fig([plotDir, 'ActiveLevelAtlas_', fileName, '.pdf']);
close

figure,
edges = -0.05:0.1:1.05;
distAct = histcounts(activeLevelAfter(activeTagBefore), edges)/sum(activeTagBefore);
distInAct = histcounts(activeLevelAfter(~activeTagBefore), edges)/sum(~activeTagBefore);
diffPdf = cumsum(distInAct) - cumsum(distAct);
[~, thresIdx] = max(diffPdf);
activeThres = edges(thresIdx + 1);

bar(edges(1:end-1)+0.05, [distInAct; distAct]', 'hist');
hold on,
plot([activeThres, activeThres], [0, 0.6], '--r');
legend({'non-active before', 'active before'}, 'Location', 'northwest');
xlabel('activation level after ablation');
ylabel('percentage');
export_fig([plotDir, 'ActiveLevelDistr_', fileName, '.pdf']);
close

movefile([tempDatDir, fileName, '.mat'], [tempDatDir, fileName, '_full', '.mat']);
save([tempDatDir, fileName, '_full', '.mat'], 'activeThres', '-append');
movefile([tempDatDir, fileName, '_tmp', '.mat'], [tempDatDir, fileName, '.mat']);

end
