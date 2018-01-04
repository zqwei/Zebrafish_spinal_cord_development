addpath('../Func');
setDir;
inputBase = '../Data';
outputBase = '../../Data';

for i = 3:24
    inputFolder = [inputBase '\' dataset{i}];
    outputFile = [outputBase '\' dataset{i} '\profile.mat'];
    if exist([inputFolder '\profile_mnx_r8.mat'], 'file')
        load([inputFolder '\profile_mnx_r8.mat'], 'profile_all');
        mnx_level_func = profile_all;
        save(outputFile, 'mnx_level_func', '-append');
    end

    if exist([inputFolder '\birthtime.mat'], 'file')
        load([inputFolder '\birthtime.mat'], 'birthtime', 'siblings');
        save(outputFile, 'birthtime', 'siblings', '-append');
    end
    
    if exist([inputFolder '\islet.mat'], 'file')
        load([inputFolder '\islet.mat'], 'islet');
        save(outputFile, 'islet', '-append');
    end
end