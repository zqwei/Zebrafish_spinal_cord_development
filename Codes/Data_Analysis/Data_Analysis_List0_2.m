% 
% Data Analysis List
% 
% -------------------------------------------------------------------------
%
% version 1.2
% preprocessing of orginal data

% function Data_Analysis_List0_2

    PADir                 = '../Preprocessing_Analysis/';
    addpath('../Func');
    % load data
    setDir;
    
    lenWindow             = 511;
    numOrder              = 9;
    
    fid                   = fopen([PADir 'fileList.dat'], 'w');
    
    for nFile             = 1:length(fileDirNames) %#ok<SUSENS>
        fileDirName       = fileDirNames{nFile};
        fileName          = fileNames{nFile}; %#ok<SUSENS>
        dirImageData      = [fileDirName '/'];
        
        load ([dirImageData 'dff.mat'], 'baseline')
        load ([dirImageData 'profile.mat'], 'profile_all')
        
        [numNeuron, numTime] = size(profile_all); %#ok<SUSENS>
        fwrite(fid, [fileName, ', ' num2str(numNeuron)])
        
        nNeuronDat        = zeros(numTime, 3);
        
        for nNeuron       = 1:numNeuron
            nNeuronDat(:, 1) = profile_all(nNeuron, :);
            nNeuronDat(:, 2) = baseline(nNeuron, :); %#ok<SUSENS>
            nNeuronDat(:, 3) = sgolayfilt(profile_all(nNeuron,:), numOrder, lenWindow);
            csvwrite([PADir fileName '_' num2str(nNeuron, '%03d')], nNeuronDat);
        end
    end
    
    fclose(fid);
% end