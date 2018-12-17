%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure materials for 2C, 2D, 2E, 2F, 2G, 2H, 2I
%
% # of active/unfactored neurons
% # of communities over time
% size of communities - scattered plot 
% activation time vs. segment
% 
% multi-fish statistics aligned by 50% of fracActNeurons
% Example fish statistics - 2C, 2D, 2E
% And a fit for all the curves - 2G, 2H, 2I
% mnx vs. activation time - 2F
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
% 

control_datasets = [3, 4, 7, 10, 12, 15, 16];

tag =  'patternTime_grpstats'; %'actTime'; %
h = figure('Position', [0, 0, 1500, 200]);
hold on;

addpath('../Func');
setDir;

mColor = lines(numel(control_datasets));
for i = 1:numel(control_datasets)
    Figure_2(h, control_datasets(i), tag, mColor(i, :));
end
setPrint(8*6, 6, [plotDir 'Figure_2G-I' tag '_hp'], 'pdf')
