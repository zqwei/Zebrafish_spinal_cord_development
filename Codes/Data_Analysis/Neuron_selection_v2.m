for nFile = 1:3
    addpath('../Func');
    setDir;
    
    fileDirName   = fileDirNames{nFile};
    fileName      = fileNames{nFile};
    
    lenWindow     = 511;
    numOrder      = 9;

    dirImageData  = [fileDirName '/'];

    load ([dirImageData 'profile.mat']) % profile_all, timepoints, nCells, tracks_smoothed, side
    load([tempDatDir, fileName, '.mat'], 'leafOrder', 'slicedIndex');

    tracks        = tracks_smoothed;
    
    timeStep      = 1200;
    numT          = size(tracks, 2);
    timeEnd       = (floor(numT/timeStep) - 2)*timeStep;    

    tracks        = tracks(slicedIndex, timeStep:timeEnd+timeStep, :);
    tracks        = tracks(leafOrder, :, :);

    save([tempDatDir, fileName, '.mat'], 'tracks', '-append');    
end