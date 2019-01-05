%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell_Segmentation - segment cell structures
% based on reference time point before ablation, generate Atlas file
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
function Cell_Segmentation(nFile)
addpath('../Func');
setDir;

backgroundThres = 130;


templateFile = recoverFilenameFromPattern([inputFolderList{nFile, 1} '/SPM00/TM??????/SPM00_TM??????_CM01_CHN00.klb'], inputFolderList{nFile, 2});
[templateFolder, ~, ~] = fileparts(templateFile);
templateImage = readImage(templateFile);
writeImage(templateImage, [templateFolder '\template.tif']);
cd('Tracking_GMM_project-v0.2.7-win64\bin\');
system(['ProcessStack.exe ' templateFile(1:end-4)  ' 2 2 ' num2str(backgroundThres) ' 74'], '-echo');
system(['ProcessStack.exe ' templateFile(1:end-4)  '_seg_conn74_rad2.bin 2 50'], '-echo');
cd('..\..\')
testSegmentation([templateFile(1:end-4) '_seg_conn74_rad2.bin_tau2.klb']);


if ~exist([inputFolderList{nFile, 1} '/Markers'], 'dir')
    mkdir([inputFolderList{nFile, 1} '/Markers']);
end
atlasImage = readImage(recoverFilenameFromPattern([inputFolderList{nFile, 1} '/SPM00/TM??????/SPM00_TM??????_CM01_CHN01.klb'], floor(inputFolderList{nFile, 2}/chn01_interval)*chn01_interval));
atlasImage = permute(atlasImage, [2, 1, 3]);
writeImage(atlasImage, [inputFolderList{nFile, 1} '/Markers/Atlas.tif']);


