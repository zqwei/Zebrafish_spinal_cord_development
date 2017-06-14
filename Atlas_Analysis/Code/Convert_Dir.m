addpath('../Func');
setDir;
inputFolder = 'D:\Code\Ziqiang_git\TempDat';
outputBase = 'D:\Code\atlas3D\Data';
for i = 1:numel(dataset)
    load([inputFolder '\' dataset{i} '.mat']);
    sum(sum(activeNeuronMat, 2)==0)
    outputFolder = [outputBase '\' dataset{i}];
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
    copyfile([inputFolder '\' dataset{i} '.mat'], [outputFolder '\data.mat']);
    copyfile([inputFolder '\EV_' dataset{i} '.mat'], [outputFolder '\EV.mat']);
    copyfile([inputFolder '\FALONO_' dataset{i} '.mat'], [outputFolder '\FALONO.mat']);
    copyfile([inputFolder '\LONOLoading_' dataset{i} '.mat'], [outputFolder '\LONOLoading.mat']);
end