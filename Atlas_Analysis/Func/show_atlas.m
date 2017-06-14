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
    xlim([0 10]);
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
    xlim([0 10]);
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
    xlim([0 10]);
    ylim([0, max(z)]);
    xlabel('AP');
    ylabel('DV');
end