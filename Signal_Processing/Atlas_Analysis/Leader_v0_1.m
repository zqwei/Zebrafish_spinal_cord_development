%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% v0_1.  Test valid RSquare threshold for halfEVTime fit
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
% 

function [numOutlier, numFit] = Leader_v0_1(nFile, RSquareThres) 
addpath('../Func');
setDir;
fileName          = fileNames{nFile};



load([tempDatDir, 'EV_', fileName, '.mat'], 'halfEVTime', 'firstActTime', 'validFitIndex', 'RSquare', 'halfEVTimeCI');

% pattern time and active time
patternTime = halfEVTime;
patternTime(RSquare<=RSquareThres | ~validFitIndex)= NaN;
patternTimeBound = halfEVTimeCI(:, 2);
patternTimeBound(RSquare<=RSquareThres | ~validFitIndex) = NaN;
activeTime = firstActTime/60;
numOutlier = sum((activeTime - patternTimeBound)>0.2);
numFit = sum(~isnan(patternTime));
end
