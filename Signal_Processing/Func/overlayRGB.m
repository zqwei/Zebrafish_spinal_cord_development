%
% overlayRGB.m
%
% Brief introduction of the usage
% 
%
% version 1.0
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%
%
%







function overlay = overlayRGB(image, colorMask)
    overlay      = colorMask;
    
    if [size(image), 3] ~= size(colorMask);
        disp('image and colorMask mismatch');
        return;
    end
    
    image = double(image);
    image = (image - min(image(:)))/(max(image(:)) - min(image(:)));
    
    if numel(size(image))==2 %2d RGBimage
        for i = 1:size(image, 1)
            for j = 1:size(image, 2)
                if (colorMask(i, j, :)==0)
                    overlay(i, j ,:) = image(i, j);
                end
            end
        end
    elseif numel(size(image))==3 %3d RGB stack
        for i = 1:size(image, 1)
            for j = 1:size(image, 2)
                for k = 1:size(image, 3)
                    if (colorMask(i, j, k, :)==0)
                        overlay(i, j ,k, :) = image(i, j, k);
                    end
                end
            end
        end
    else
        disp('image type not supported');
    end

