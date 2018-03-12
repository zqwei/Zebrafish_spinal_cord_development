%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 1D plot 3 example long-term traces
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
% 
addpath('../Func')
setDir;
load([fileDirNames{3} '/profile.mat'], 'timepoints');
load([tempDatDir '/' fileNames{3} '.mat'], 'dff');

cell_id = [37, 16, 70];
% cell_id = [13, 14, 78, 79, 90];

timepoints = (1:size(dff, 2)) * 0.25/3600;

figure('Position', [0, 100, 1000, 300]);
plot(timepoints, dff(cell_id, :)+ repmat(-(0:(numel(cell_id)-1))'*2, 1, numel(timepoints)));
xlim([0, max(timepoints)]);
colormap(lines)
set(gca, 'yticklabel', []);
set(gca, 'ytick', []);
set(gca, 'ycolor', 'w');
box off
setPrint(20, 6, [plotDir 'Figure_1D Example Traces'] , 'pdf');
close

