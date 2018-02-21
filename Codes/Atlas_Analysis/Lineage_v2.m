%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lineage analysis 2: segmental analysis
%
% birth order and activation order within segment
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%

function Lineage_v2(nFile)
addpath('../Func');
setDir;

fileName          = fileNames{nFile};
load([tempDatDir, 'Leader_', fileName, '.mat'], 'activeTime', 'birthtime');
load([tempDatDir, fileName, '.mat'], 'timePoints', 'new_x', 'new_y', 'mnx');
actOrder = nan(numel(new_x), 1);
birthOrder = nan(numel(new_x), 1);

x = new_x;
y = new_y;
for seg = floor(min(x-0.5)):ceil(max(x-0.5))
    currentSegLeft = find(y<0 & x>=seg+0.5 & x<seg+1.5 & mnx==1 & ~isnan(activeTime) & ~isnan(birthtime));
    [~, ~, ic] = unique(activeTime(currentSegLeft));
    actOrder(currentSegLeft) = ic;
    [~, ~, ic] = unique(birthtime(currentSegLeft));
    birthOrder(currentSegLeft) = ic;
    currentSegRight = find(y>0 & x>=seg+0.5 & x<seg+1.5 & mnx==1 & ~isnan(activeTime) & ~isnan(birthtime));
    [~, ~, ic] = unique(activeTime(currentSegRight));
    actOrder(currentSegRight) = ic;
    [~, ~, ic] = unique(birthtime(currentSegRight));
    birthOrder(currentSegRight) = ic;
end

count = hist3([actOrder, birthOrder], [max(actOrder), max(birthOrder)]);
figure, imagesc(count, [0 max(count(:))+1]);
set(gca, 'XTick', 1:max(birthOrder), 'YTick', 1:max(actOrder), 'YDir', 'normal');
colormap(hot);
colorbar
xlabel('birth order')
ylabel('EV order')
