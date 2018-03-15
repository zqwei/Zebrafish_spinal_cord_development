%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3B
%
% summary composition of initial pairs
% 
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org

addpath('../Atlas_Analysis/')
control_datasets = [3, 4, 7, 12, 10, 15, 16];

typeList = [];
for nFile = control_datasets
    [histPairType, actLevelBeforePattern] = Leader_v4(nFile);
    typeList = [typeList; histPairType'];
end


composition = typeList./repmat(sum(typeList, 2), 1, 3);
figure,  boxplot(composition)
set(gca, 'XTickLabel', {'L-L', 'L-N', 'N-N'});
pLN_LL = signrank(composition(:, 2), composition(:, 1), 'tail', 'right');
pLN_NN = signrank(composition(:, 2), composition(:, 3), 'tail', 'right');