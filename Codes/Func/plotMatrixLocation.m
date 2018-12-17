function newIm = plotMatrixLocation(matrix, regionSize, halfFish)
% image size is drawn proportional to regionSize
p{1} = imread('pattern0_1.tif'); % type 0
p{2} = imread('pattern1_1.tif'); % type 1
p{3} = imread('pattern2_1.tif'); % type 2

lx = size(p{1}, 1);
ly = size(p{1}, 2);
newIm = zeros(size(matrix, 1)* lx, size(matrix, 2)*ly, 3, 'uint8');

% relativeP = bsxfun(@rdivide, regionSize, sum(regionSize));
rs = floor(bsxfun(@rdivide, regionSize, sum(regionSize))*ly*3);
for nFish = 1:size(matrix, 2)
    xOffset = 0;
    for nRegion = 1:size(matrix, 1)
        newIm((1:rs(nRegion, nFish)) + xOffset, (1:ly) + (nFish-1)*ly, :) = imresize(p{matrix(nRegion, nFish)+1}, [rs(nRegion, nFish), ly]);
        xOffset = xOffset + rs(nRegion, nFish);
    end
end

if (halfFish)
    newIm = imresize(newIm, [size(matrix, 1)* lx, size(matrix, 2)*ly/2]);
end

figure, imshow(newIm);