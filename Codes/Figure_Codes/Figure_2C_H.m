%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure materials for 2C, 2D, 2E, 2F, 2G, 2H, 2I
%
% # of active/unfactored neurons
% # of communities over time
% size of communities - scattered plot 
% activation time vs. segment
% 
% Example fish statistics - 2C, 2D, 2E
% And a fit for all the curves - 2G, 2H, 2I
% mnx vs. activation time - 2F
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
% 

control_datasets = [3, 4, 7, 12, 10, 11, 15, 16];

setDir;    
h = figure('Position', [0, 0, 1500, 200]);
hold on;
for nFile = control_datasets
    Figure_2(h, nFile);
end
setPrint(8*6, 6, [plotDir 'Figure_2G-I_alignByFracAct50'], 'pdf')
