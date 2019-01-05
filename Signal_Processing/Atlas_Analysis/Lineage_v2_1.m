%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lineage analysis 2.1: segmental analysis, summarize a list of datasets
%
% birth order and activation order within segment
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%

function Lineage_v2_1(datasets)
addpath('../Func')
setDir;
birthOrderList = [];
actOrderList = [];
nSegAll = 0;
for nFish = 1:numel(datasets)
    nFile = datasets(nFish);
    [birthOrder, actOrder, nSeg] = Lineage_v2(nFile);
    birthOrderList = [birthOrderList; birthOrder];
    actOrderList = [actOrderList; actOrder];
    nSegAll = nSegAll + nSeg;
end

actOrderList(isnan(actOrderList)) = [];
birthOrderList(isnan(birthOrderList)) = [];

[r, p] = corr(actOrderList, birthOrderList);

count = hist3([actOrderList, birthOrderList], [max(actOrderList), max(birthOrderList)]);
figure, imagesc(count, [0 max(count(:))+1]);
set(gca, 'XTick', 1:max(birthOrder), 'YTick', 1:max(actOrder), 'YDir', 'normal');
colormap(hot);
colorbar
xlabel('birth order')
ylabel('activation order')
title(['#segs=' num2str(nSegAll) ', #neurons=' num2str(numel(actOrderList)) ', r=' num2str(r, '%.3f'), ', p=' num2str(p, '%.3f')]);
setPrint(8, 6, [plotDir,  'SegmentalBirthActOrder_all'], 'pdf');
