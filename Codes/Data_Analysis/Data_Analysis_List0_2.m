% 
% Data Analysis List
% 
% -------------------------------------------------------------------------
%
% version 1.2
% preprocessing of orginal data

function Data_Analysis_List0_2

    PADir                 = '../Preprocessing_Analysis/Dat/';
    addpath('../Func');
    % load data
    setDir;
    
    lenWindow             = 511;
    numOrder              = 9;
    numNeurons            = zeros(length(fileDirNames), 1);
    
    for nFile             = 1:length(fileDirNames) 
        fileDirName       = fileDirNames{nFile};
        fileName          = fileNames{nFile}; 
        dirImageData      = [fileDirName '/'];
        
        load ([dirImageData 'dff.mat'], 'baseline')
        load ([dirImageData 'profile.mat'], 'profile_all')
        
        [numNeuron, numTime] = size(profile_all); 
        numNeurons(nFile)    = numNeuron;
        
        nNeuronDat        = zeros(numTime, 3);
        
        for nNeuron       = 1:numNeuron
            nNeuronDat(:, 1) = profile_all(nNeuron, :);
            nNeuronDat(:, 2) = baseline(nNeuron, :); 
            nNeuronDat(:, 3) = sgolayfilt(profile_all(nNeuron,:), numOrder, lenWindow);
            save([PADir fileName '_' num2str(nNeuron, '%03d') '.mat'], 'nNeuronDat');
        end
    end
    
    save([PADir 'fileList.mat'], 'fileNames', 'numNeurons')
end