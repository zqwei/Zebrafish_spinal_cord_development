addpath('../Func');
setDir;
inputFolder = '../../TempDat';
outputBase = '../Data';

exMetrFolder = 'D:\Code\atlas3D\Data';
for i = 23:24
    load([inputFolder '\' dataset{i} '.mat']);
    % number of never-active neurons
    % sum(sum(activeNeuronMat, 2)==0)
    outputFolder = [outputBase '\' dataset{i}];
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
    copyfile([inputFolder '\' dataset{i} '.mat'], [outputFolder '\data.mat']);
    copyfile([inputFolder '\EV_' dataset{i} '.mat'], [outputFolder '\EV.mat']);
    copyfile([inputFolder '\FALONO_' dataset{i} '.mat'], [outputFolder '\FALONO.mat']);
    copyfile([inputFolder '\LONOLoading_' dataset{i} '.mat'], [outputFolder '\LONOLoading.mat']);
    
    if exist([inputFolder '\LONOLoading_v_0_1' dataset{i} '.mat'], 'file')
            copyfile([inputFolder '\LONOLoading_v_0_1' dataset{i} '.mat'], [outputFolder '\LONOLoading_v_0_1.mat']);
    end
    
    exfolder = [exMetrFolder '\' dataset{i}];
    if exist([exfolder '\islet.mat'], 'file')
        copyfile([exfolder '\islet.mat'], [outputFolder '\islet.mat']);
    end
    if exist([exfolder '\profile_mnx_r8.mat'], 'file')
        copyfile([exfolder '\profile_mnx_r8.mat'], [outputFolder '\profile_mnx_r8.mat']);
    end
    if exist([exfolder '\birthtime.mat'], 'file')
        copyfile([exfolder '\birthtime.mat'], [outputFolder '\birthtime.mat']);
    end
end