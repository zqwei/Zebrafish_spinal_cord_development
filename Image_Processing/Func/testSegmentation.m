%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% testSegmentation.m
% output segmentation result as mdf file
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%
function testSegmentation(inputStackName)
[inputFolder, ~, ~] = fileparts(inputStackName);
stack = readImage(inputStackName);
centers = regionprops(stack, 'Centroid');

pointList = cell2mat(struct2cell(centers)');
pointList(:, 3) = round(pointList(:, 3));
writeMDF([inputFolder '/segmentation.mdf'], pointList);
