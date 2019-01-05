%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre-processing Step 2 -- Image segmentation
% Segmentation performed by TGMM cell tracking, see
%        Amat, Fernando, et al. "Fast, accurate reconstruction of cell lineages from large-scale fluorescence microscopy data." Nature methods 11.9 (2014): 951.
% Ouput: mdf file for visluzation by MTrackJ:
%        https://imagescience.org/meijering/software/mtrackj/
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
%


function Pre_Processing_s2(inputFolder, backgroundThres)
addpath('../Func');
load([inputFolder 'signal\periods.mat'], 'nPeriod', 'periodStart');

templateFile = [inputFolder 'signal\Period_TM' num2str(nPeriod-2, '%.6d') '.klb'];
[templateFolder, ~, ~] = fileparts(templateFile);
templateImage = readImage(templateFile);
writeImage(templateImage, [templateFolder '\template.tif']);
cd('Tracking_GMM_project-v0.2.7-win64\bin\');
system(['ProcessStack.exe ' templateFile(1:end-4)  ' 2 2 ' num2str(backgroundThres) ' 74'], '-echo');
system(['ProcessStack.exe ' templateFile(1:end-4)  '_seg_conn74_rad2.bin 2 50'], '-echo');
cd('..\..\')

inputStackName = [templateFile(1:end-4) '_seg_conn74_rad2.bin_tau2.klb'];
testSegmentation(inputStackName);
stack = readImage(templateFile);
centers = regionprops(stack, 'Centroid');

pointList = cell2mat(struct2cell(centers)');
pointList(:, 3) = round(pointList(:, 3));
writeMDF([inputStackName '_truth.mdf'], pointList);
