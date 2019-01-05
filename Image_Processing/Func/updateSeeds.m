%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update location of the seed point
% ratio3D: set to [] for 2D tracking, otherwise use floor(windowSize/ratio)
% as window size in the axial direction
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%
function [results, label] = updateSeeds(fileName, seeds, oldLabel, windowSize, ratio3D, maxDistance)


image = readImage(fileName);
for i = 1:size(image, 3)
    image(:, :, i) = medfilt2(squeeze(image(:, :, i)));
end
mask = zeros(size(image));

results = seeds;
label = oldLabel;

nCells = size(seeds, 1);

label_merge = false(nCells, 1);

%% region growing based recenter
for i = 1:nCells
    if label(i) == 0
        continue;
    end
    cx = round(seeds(i, 2));
    cy = round(seeds(i, 1));
    cz = round(seeds(i, 3));
    wx = min([windowSize, cx-1, size(image, 1)-cx]);
    wy = min([windowSize, cy-1, size(image, 2)-cy]);
    if (wx < windowSize || wy < windowSize)
        label(i) = 0;
        disp(['    cell ' num2str(i) ' eliminated, too close to xy boundary']);
        continue;
    end
    if ~isempty(ratio3D)
        wz = min([floor(windowSize/ratio3D), cz-1, size(image, 3)-cz]);
        if wz < 1
            label(i) = 0;
            disp(['    cell ' num2str(i) ' eliminated, too close to z boundary']);
            continue;
        end
    else
        wz = 0;
    end
    
    localWindow = image(cx-wx:cx+wx, cy-wy:cy+wy, cz-wz:cz+wz);

    [local_cx, local_cy, local_cz, local_mask] = recenterRG(wx+1, wy+1, wz+1, localWindow, ratio3D);
    
    results(i, 2) = cx - wx -1 + local_cx;
    results(i, 1) = cy - wy -1 + local_cy;
    if ~isempty(ratio3D)
        results(i, 3) = cz - wz -1 + local_cz;
    end
    
    
    %% section for debugging
%     subplot(2, 1, 1),
%     imagesc(local_mask(:,:, wz+1));
%     hold on
%     scatter(wy+1, wx+1, '*');
%     scatter(local_cy, local_cx, 'o');
%     title(['cell ' num2str(i)]);
%     hold off
%     subplot(2, 1, 2),
%     imagesc(localWindow(:,:, wz+1));
%     hold on
%     scatter(wy+1, wx+1, '*');
%     scatter(local_cy, local_cx, 'o');
%     title(['cell ' num2str(i)]);
%     hold off
    
    
    %% split if overlap in mask is detected
    [mx, my, mz] = ind2sub(size(local_mask), find(local_mask(:)==1));
    global_ind = sub2ind(size(mask), cx - wx -1 + mx, cy - wy -1 + my, cz - wz -1 + mz);
    if any(mask(global_ind)>0)
        cells_to_merge =  unique(mask(global_ind));
        cells_to_merge(cells_to_merge==0) = [];
        for k = 1:numel(cells_to_merge)
            disp(['detected merging of cell ' num2str(i) ' with cell ' num2str(cells_to_merge(k))]);
            mask(global_ind) = i;
            [centers_new, mask, valid_split] = recenter_clustering(mask, seeds, image, i, cells_to_merge(k), ratio3D);
            results(i, :) = centers_new(1, :);
            results(cells_to_merge(k), :) = centers_new(2, :);
            if (valid_split)
                label_merge(i) = 1;
                label_merge(cells_to_merge(k)) = 1;
            else
                vec_i = results(i, 1:2) - seeds(i, 1:2);
                vec_merge = results(cells_to_merge(k), 1:2) - seeds(cells_to_merge(k), 1:2);
                if (norm(vec_i) > norm(vec_merge))
                    disp(['    invalid merge detected, kill cell ' num2str(i) ]);
                    label(i) = 0;
                    results(cells_to_merge(k), :) = seeds(cells_to_merge(k), :);
                else
                    disp(['    invalid merge detected, kill cell ' num2str(cells_to_merge(k)) ]);
                    label(cells_to_merge(k)) = 0;
                    results(i, :) = seeds(i, :);
                end
            end
        end
        
    else
        mask(global_ind) = i;
    end


    
end


%% sanity check
% if after updating and merging, the cells still move more than
% maxDistance, keep the original location

for i = 1:nCells
    if label(i)==0
        continue;
    end
    vec = seeds(i, :) - results(i, :);
    dist = sqrt(vec(1)^2 + vec(2)^2);
    if (dist>maxDistance || abs(vec(3)) > 1)
        if ~label_merge(i)
            disp(['    cell ' num2str(i) ' warning: moved too far across one frame']);
%             results(i, :) = seeds(i, :);
        end
    end
end
