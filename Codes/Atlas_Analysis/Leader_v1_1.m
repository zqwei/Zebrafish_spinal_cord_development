%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.1. Plot metrics on atlas
%
% all metrics prestored in "listLeaderMetrics"
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
% 

function Leader_v1_1(nFile) 
addpath('../Func');
setDir;
fileName          = fileNames{nFile};

% calculate factor size
load([tempDatDir, fileName, '.mat'], 'new_x', 'new_y', 'new_z', 'side', 'timePoints', 'mnx');
load([tempDatDir, 'Leader_' fileName, '.mat']);
load([tempDatDir, 'EV_', fileName, '.mat'], 'halfActTime', 'halfEVTime');


% calculate other metrics
diffHalfTime = halfEVTime - halfActTime;

% plot atlas with different metrics
x = new_x;
y = new_y;
z = new_z;
for it = 1:numel(listLeaderMetrics)
    me = eval(listLeaderMetrics{it});
    show_atlas(x, y, z, 1:max(x), me, mnx>0);
    export_fig([plotDir '/' listLeaderMetrics{it} '_' fileName '.pdf']);
    close
end
end

function show_atlas(x, y, z, ra, profile, label)
% plot the atlas with coordinates [x, y, z, ra]
% with color showing profile
% label=1 shown as filled, label=0 shown as +
x(isnan(profile)) = [];
y(isnan(profile)) = [];
z(isnan(profile)) = [];
label(isnan(profile)) = [];
profile(isnan(profile)) = [];

figure,
    % dorsal view
    subplot(3, 1, 1)
    hold on
    scatter(x(label), y(label), [], profile(label), 'filled');
    scatter(x(~label), y(~label), [], profile(~label), '+');
    for i = 0:numel(ra)+1
        plot([i, i], [min(y), max(y)], '--k');
    end
    caxis([min(profile), max(profile)]);
    colorbar
    hold off
    title('dorsal view');
    xlabel('AP');
    ylabel('ML');
    xlim([0 numel(ra)+1]);
    ylim([min(y), max(y)]);
    
    % side view - left
    subplot(3, 1, 2)
    hold on
    scatter(x(label & y<0), z(label & y<0), [], profile(label & y<0), 'filled');
    scatter(x(~label & y<0), z(~label & y<0), [], profile(~label & y<0), '+');
    for i = 0:numel(ra)+1
        plot([i, i], [0, max(z)], '--k');
    end
    caxis([min(profile), max(profile)]);
    colorbar
    hold off
    title('side view - right');
    xlim([0 numel(ra)+1]);
    ylim([0, max(z)]);
    xlabel('AP');
    ylabel('DV');
    
    
    % side view - right
    subplot(3, 1, 3)
    hold on
    scatter(x(label & y>0), z(label & y>0), [], profile(label & y>0), 'filled');
    scatter(x(~label & y>0), z(~label & y>0), [], profile(~label & y>0), '+');
    for i = 0:numel(ra)+1
        plot([i, i], [0, max(z)], '--k');
    end
    caxis([min(profile), max(profile)]);
    colorbar
    hold off
    title('side view - left');
    xlim([0 numel(ra)+1]);
    ylim([0, max(z)]);
    xlabel('AP');
    ylabel('DV');
end