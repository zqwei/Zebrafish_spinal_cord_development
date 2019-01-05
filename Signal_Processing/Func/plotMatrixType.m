function newIm = plotMatrixType(matrix, halfFish)
p{1} = imread('pattern0_1.tif'); % type 0
p{2} = imread('pattern1_1.tif'); % type 1
p{3} = imread('pattern2_1.tif'); % type 2

lx = size(p{1}, 1);
ly = size(p{1}, 2);
newIm = zeros(size(matrix, 1)* lx, size(matrix, 2)*ly, 3, 'uint8');
for i = 1:size(matrix, 1)
    for j = 1:size(matrix, 2)
        newIm((1:lx)+(i-1)*lx, (1:ly) + (j-1)*ly, :) = p{matrix(i, j)+1};
    end
end
if (halfFish)
    newIm = imresize(newIm, [size(matrix, 1)* lx, size(matrix, 2)*ly/2]);
end

figure, imshow(newIm);