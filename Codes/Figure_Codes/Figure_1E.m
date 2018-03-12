%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 1E plot 3 example time windows
%
% All traces in grey, example traces highlighted in color
% window length: 1 min
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
% 

addpath('../Func')
setDir;
load([tempDatDir '/' fileNames{3} '.mat'], 'dff', 'new_x', 'new_y', 'timePoints');

example_windows = [10, 120, 191];
cell_id = [37, 16, 70];

color_table = lines(numel(cell_id));
spacing = 0.3;
length_time = 240; % in time points
[~, neworder] = sort(new_x);
neworder = [neworder(new_y(neworder)<0); neworder(new_y(neworder)>0)];
boundary = -sum(new_y<0) * spacing-spacing/2;
for i = 1:3  
    current_window = timePoints(example_windows(i)):timePoints(example_windows(i))+length_time;
    figure('position', [0, 100, 300, 500]),
    hold on
    dff_plot = dff(:, current_window);
%     % if using zscore
%     dff_plot = zscore(dff_plot, [], 2);
    
    plot((1:numel(current_window))/4, dff_plot(neworder, :) + repmat(-(1:size(dff, 1))'*spacing, 1, numel(current_window)), 'color', [.5, .5, .5]);
    xlim([0, length_time * 0.25]);
    plot([0, length_time * 0.25], [boundary, boundary], '--k');
    for j = 1:numel(cell_id)
        plot((1:numel(current_window))/4, dff_plot(cell_id(j), :)+ repmat(-find(neworder==cell_id(j))*spacing, 1, numel(current_window)), 'color', color_table(j, :), 'linewidth', 1.5);
    end
    box on
    set(gca, 'ytick', []);
    set(gca, 'yticklabel', []);
%     set(gca, 'xtick', 0:5);
    ylim([floor(-numel(neworder)*spacing), spacing]);
    setPrint(6, 10, [plotDir 'Figure_1E Example Window_' num2str(example_windows(i))] , 'pdf');
    close;
end