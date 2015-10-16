%
% overlayRGB.m
%
% Brief introduction of the usage
% 
% impose color mask at the coordinates accordinates with the colors in
% colorTable.
% Input:
% inputFileName:    the original image to impose on
% coordinates:      N x 3 matrix specifying the x,y,z coordinates of cells
% colorTable:       N x 3 matrix specifying the r,g,b color of cells
%
% version 1.0
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%
%
%

function exportMaskRadius(inputFileName, outputFilePath, coordinates, colorTable, radiusTable)

    ratio                = 6; % what is radius and ratio
    jpegWorkers          = 12;

    [~, ~, fileType]     = fileparts(inputFileName);
    if strcmp(fileType, '.jp2')
        overlay          = readJStack(inputFileName, jpegWorkers);
    elseif strcmp(fileType, '.tif') || strcmp(fileType, '.tiff')
        overlay          = readTIFF(inputFileName);
    end

    stackDimension       = size(overlay);
    overlay_xyProjection = squeeze(max(overlay, [], 3));
    overlay_yzProjection = squeeze(max(overlay, [], 1));
    overlay_xzProjection = squeeze(max(overlay, [], 2));

    overlay              = double(overlay);
    overlay              = repmat((overlay-min(overlay(:)))/(max(overlay(:))-min(overlay(:))), [1 1 1 3]);
    coordinates          = round(coordinates);

    rgbStack = zeros([stackDimension, 3]);

    for i = 1:size(coordinates, 1)
        radius       = ceil(radiusTable(i)*20);
        [rr, cc, zz] = meshgrid(-radius:radius, -radius:radius, round(-radius/ratio):round(radius/ratio));
        [cmx, cmy, cmz] = ind2sub([2*radius+1, 2*radius+1, 2*round(radius/ratio)+1], find(rr.^2 + cc.^2 + (ratio*zz).^2<= radius^2));
        for j = 1:length(cmx)
            overlay(cmx(j) - radius - 1 + coordinates(i, 2) , cmy(j) - radius -1 + coordinates(i, 1), cmz(j) - round(radius/ratio) - 1 + coordinates(i, 3), :) = colorTable(i, :);
            rgbStack(cmx(j) - radius - 1 + coordinates(i, 2) , cmy(j) - radius -1 + coordinates(i, 1), cmz(j) - round(radius/ratio) - 1 + coordinates(i, 3), :) = colorTable(i, :);
        end
    end

    if ~exist(outputFilePath, 'dir')
        mkdir(outputFilePath);
    end

    outputFileName = [outputFilePath '/mask.tif'];

    imwrite(squeeze(rgbStack(:, :, 1, :)), outputFileName, 'tif', 'Compression', 'none');
    for i = 2:stackDimension(3)
        imwrite(squeeze(rgbStack(:, :, i, :)), outputFileName, 'tif', 'WriteMode', 'append', 'Compression', 'none');
    end
    rgbStack_xyProjection = squeeze(max(rgbStack, [], 3));
    rgbStack_yzProjection = squeeze(max(rgbStack, [], 1));
    rgbStack_xzProjection = squeeze(max(rgbStack, [], 2));
    imwrite(rgbStack_xyProjection, [outputFileName '_xyProjection.tif'], 'tif', 'Compression', 'none');
    imwrite(imresize(rgbStack_yzProjection, [stackDimension(2), stackDimension(3)*ratio]), [outputFileName '_yzProjection.tif'], 'tif', 'Compression', 'none');
    imwrite(imresize(rgbStack_xzProjection, [stackDimension(1), stackDimension(3)*ratio]), [outputFileName '_xzProjection.tif'], 'tif', 'Compression', 'none');



    overlay_xyProjection = overlayRGB(overlay_xyProjection, rgbStack_xyProjection);
    overlay_yzProjection = overlayRGB(overlay_yzProjection, rgbStack_yzProjection);
    overlay_xzProjection = overlayRGB(overlay_xzProjection, rgbStack_xzProjection);

    overlayName = [outputFileName '_overlay.tif'];
    imwrite(squeeze(overlay(:, :, 1, :)), overlayName, 'tif', 'Compression', 'none');
    for i = 2:stackDimension(3)
        imwrite(squeeze(overlay(:, :, i, :)), overlayName, 'tif', 'WriteMode', 'append', 'Compression', 'none');
    end
    imwrite(overlay_xyProjection, [overlayName '_xyProjection.tif'], 'tif', 'Compression', 'none');
    imwrite(imresize(overlay_yzProjection, [stackDimension(2), stackDimension(3)*ratio]), [overlayName '_yzProjection.tif'], 'tif', 'Compression', 'none');
    imwrite(imresize(overlay_xzProjection, [stackDimension(1), stackDimension(3)*ratio]), [overlayName '_xzProjection.tif'], 'tif', 'Compression', 'none');

end
    