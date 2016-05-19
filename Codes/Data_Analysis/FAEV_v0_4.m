%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.  Explained variance: Factor analysis -- corrected loading matrix --
% sigmoid fit
% 
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 


function FAEV_v0_4(nFile)            
    addpath('../Func');
    setDir;    
    fileName          = fileNames{nFile}; %#ok<USENS>    
    load([tempDatDir, fileName, '.mat'], 'mnx', 'tracks', 'side', 'leafOrder', 'slicedIndex'); 
    load([tempDatDir, 'EV_' fileName, '.mat'], 'EVMat', 'halfEVTime', 'RSquare', 'halfActTime', 'validFitIndex')
    
    meantracks        = squeeze(mean(tracks, 2)); 
    RSquareThres = 0.6;    
    xTrack       = meantracks(:, 1);
    yTrack       = meantracks(:, 2);
    zTrack       = meantracks(:, 3);
    
    locationTrack = [xTrack(RSquare>=RSquareThres & validFitIndex), yTrack(RSquare>=RSquareThres & validFitIndex), zTrack(RSquare>=RSquareThres & validFitIndex)];
    halfTimeSet   = halfEVTime(RSquare>=RSquareThres & validFitIndex);
    
    [~, ~, stats]           = glmfit(locationTrack, log(halfTimeSet));
    saveIndex               = find(stats.p < 0.05);
    saveIndex(saveIndex==1) = [];
    saveIndex               = saveIndex - 1;
    [beta, ~, stats]        = glmfit(locationTrack(:, saveIndex), halfTimeSet);
    disp('linear model')
    disp(saveIndex)
    disp(beta)
    disp(stats.p)
    
    [~, ~, stats]           = glmfit(locationTrack, halfTimeSet.^2);
    saveIndex               = find(stats.p < 0.05);
    saveIndex(saveIndex==1) = [];
    saveIndex               = saveIndex - 1;
    [beta, ~, stats]        = glmfit(locationTrack(:, saveIndex), halfTimeSet.^2);
    disp('square model')
    disp(saveIndex)
    disp(beta)
    disp(stats.p)

    for nAxis = 1:3
        [rho, pval] = corr(locationTrack(:, nAxis), halfTimeSet, 'type', 'Spearman')
    end
    
    disp('-----------------')
    
end