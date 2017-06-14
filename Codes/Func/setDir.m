%%%%% Directory variables

dataFileReadMe;

fileNames          = cell(length(fileDirNames),1);

for nName          = 1:length(fileDirNames)
    [~, fileNames{nName}, ~]   = fileparts(fileDirNames{nName});
end
%%%%% Add command for warning off (especially for file loads)
warning('off','all')


plotDir             = '../../Plots/';
plotNetDir          = '../../NetworkDynamicsPlots/';

tempDatDir          = '../../TempDat/';
<<<<<<< HEAD
tempDatNetDir       = '../../NetworkDynamicsTempDat/';
=======

%if ~exist(tempDatDir, 'dir')
%    mkdir(tempDatDir)
%end
>>>>>>> master
