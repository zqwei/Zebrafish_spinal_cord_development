%
% readTIFF.m
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




function image = readTIFF(file, firstSlice, lastSlice)
 
    info = imfinfo(file);

    if nargin < 3; lastSlice = numel(info); end    
    if nargin < 2; firstSlice = 1; end
    
    numImages = 1 + lastSlice - firstSlice;
    
    bpp = info(1).BitDepth;
    
    if bpp == 1
        dtype = 'logical';
    elseif bpp <= 8
        dtype = 'uint8';
    elseif bpp <= 16
        dtype = 'uint16';
    elseif bpp <= 32
        dtype = 'single';
    else % use double
        dtype = 'double';
    end
    
    image = zeros(info(1).Height, info(1).Width, numImages, dtype);
    
    for i = firstSlice:lastSlice; image(:,:, 1 + i - firstSlice) = imread(file, i); end
end
