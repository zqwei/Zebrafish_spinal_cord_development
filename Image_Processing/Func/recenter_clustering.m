%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Split overlapping detected objects
% use k-means clustering to split overlapping
% if cell center moved too far across time, mark as invalid
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%
function [centers_new, mask, valid] = recenter_clustering(mask, seeds, image, i, j, ratio)



minNeightbor = 10; % theshold for maximum movement of cell
centers = seeds([i; j], :);
centers_new = centers;

valid = 1;

if (round(centers(1, 3))==round(centers(2, 3))) || isempty(ratio)
    cz = round(centers(1, 3));
    [cx, cy] = ind2sub(size(mask(:, :, 1)), find(mask(:, :, cz)==i | mask(:, :, cz)==j));
    [idx, centers_new(:, 1:2)] = kmeans([cy, cx], 2, 'start', centers(:, 1:2), 'emptyaction', 'drop');
    if any(isnan(centers_new(:))) %empty clusters
        valid = 0;
        disp('empty cluster detected');
        return;
    end
    for k = 1:numel(idx)
        if idx(k)==1
            mask(cx(k), cy(k), cz) = i;
        else
            mask(cx(k), cy(k), cz) = j;
        end
    end

else
    [cx, cy, cz] = ind2sub(size(mask), find(mask==i | mask==j));
    centers(:, 3) = (centers(:, 3) - 1) * ratio + 1;
    cz = (cz - 1) * ratio + 1;
    [idx, centers_new] = kmeans([cy, cx, cz], 2, 'start', centers, 'emptyaction', 'drop');
    if any(isnan(centers_new)) % case of empty clusters
        valid = 0;
        disp('empty cluster detected');
        return;
    end
    centers_new(:, 3) = (centers_new(:, 3) - 1)/ratio + 1;
    cz = (cz - 1)/ratio +1;
    for k = 1:numel(idx)
        if idx(k)==1
            mask(cx(k), cy(k), round(cz(k))) = i;
        else
            mask(cx(k), cy(k), round(cz(k))) = j;
        end
    end
end


vec = centers_new(1, :) - centers_new(2, :);
if isempty(ratio)
    dist = sqrt(vec(1)^2 + vec(2)^2);
else
    dist = sqrt(vec(1)^2 + vec(2)^2 + (vec(3)*ratio)^2);
end
if (dist<minNeightbor)
    valid = 0;
%     imagesc(image(min(cx):max(cx), min(cy):max(cy), round(mean(cz))));
%     hold on
%     scatter(centers(:, 1) - min(cy), centers(:, 2) - min(cx), '*');
%     scatter(centers_new(:, 1) - min(cy), centers_new(:, 2) - min(cx), 'o');
%     hold off
%     imagesc(mask(min(cx):max(cx), min(cy):max(cy), round(mean(cz))));
%     hold on
%     scatter(centers(:, 1) - min(cy), centers(:, 2) - min(cx), '*');
%     scatter(centers_new(:, 1) - min(cy), centers_new(:, 2) - min(cx), 'o');
%     hold off
end