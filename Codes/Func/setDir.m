%%%%% Directory variables

dataFileReadMe;

fileNames          = cell(length(fileDirNames),1);

for nName          = 1:length(fileDirNames)
    [~, fileNames{nName}, ~]   = fileparts(fileDirNames{nName});
end
%%%%% Add command for warning off (especially for file loads)
warning('off','all')

%%%%% Add Path for Plots
% addpath('../../../../Empirical_Data_Analysis_Code/Plots/');
% addpath('../../../Empirical_Data_Analysis_Code/GraphViz/');

%%%%% Add Path for LDSI Fit
addpath('../../../../Empirical_Data_Analysis_Code/LDSI/Release_LDSI_v3/');
% addpath('../../../../Empirical_Data_Analysis_Code/x3d_version3g/');



% dirImageData      = [fileDirName, '/Data/'];
% TempDat           = [fileDirName, '/Fitted_result/'];

plotDir             = '../../Plots/';

if ~exist(plotDir, 'dir')
    mkdir(plotDir)
end

tempDatDir          = '../../TempDat/';

if ~exist(tempDatDir, 'dir')
    mkdir(tempDatDir)
end