%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recenter image based on region growing method
% Background segmentation using otsu threshold
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%
function [cx_new, cy_new, cz_new, mask] = recenterRG(cx, cy, cz, image, ratio)

thresF = 0.85;
image = double(image);
image = (image - min(image(:))+1)/(max(image(:)) - min(image(:))+1);
mask = zeros(size(image));
% mask = 0: not explored
% mask = 1: foreground
% mask = 2: to be explored
% mask = -1: background
mask(cx, cy, cz) = 2;
iter = 0;
if isempty(ratio)
    ratio = 0;
end
while sum(sum(sum(mask==2))) > 0
    iter = iter + 1;
    if iter > numel(mask)
        disp('stuck in dead end');
    end
    candidatePool = find(mask==2);
    [px, py, pz] = ind2sub(size(mask), candidatePool);
    dist_pool = (px-cx).^2 + (py-cy).^2 + ((pz-cz)*ratio).^2;
    [~, idx] = min(dist_pool);
    currentCandidate = candidatePool(idx);
    if image(currentCandidate) > max(0, mean(image(mask==1)) * thresF)
        mask(currentCandidate) = 1;
        [tx, ty, tz] = ind2sub(size(image), currentCandidate);
        if tx > 1 && mask(tx-1, ty, tz) == 0
            mask( tx-1, ty, tz) = 2;
        end
        if tx < size(image, 1) && mask(tx+1, ty, tz) == 0
            mask(tx+1, ty, tz) = 2;
        end
        if ty > 1 && mask(tx, ty-1, tz) == 0
            mask(tx, ty-1, tz) = 2;
        end
        if ty < size(image, 2) && mask(tx, ty+1, tz) == 0
            mask(tx, ty+1, tz) = 2;
        end
        if tz > 1 && mask(tx, ty, tz-1) == 0
            mask(tx, ty, tz-1) = 2;
        end
        if tz < size(image, 3) && mask(tx, ty, tz+1) == 0
            mask(tx, ty, tz+1) = 2;
        end
    else
        mask(currentCandidate) = -1;
    end
end
mask(mask~=1) = 0;
result = regionprops(mask, image, 'Centroid');
cx_new = result.Centroid(2);
cy_new = result.Centroid(1);
if numel(result.Centroid)==3
    cz_new = result.Centroid(3);
else
    cz_new = [];
end

