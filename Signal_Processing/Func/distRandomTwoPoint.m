% 
% Get radius of a cluster of data (shuffled data)
%
% Input:
% 
% X        --- N x m
% numStim  --- number of stimulation
% 
% -------------------------------------------------------------------------
% version 1.0
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 

function meanDist  = distRandomTwoPoint(X, numStim)

    numNeuron      = size(X,1);
    randMat        = ceil(rand(numStim, 2) * numNeuron);
    allDist        = arrayfun(@(tIndex) getRadius(X(randMat(tIndex,:),:)),...
                                1:numStim, 'UniformOutput', false);
    meanDist       = mean(cell2mat(allDist));
    
    % hist(cell2mat(allDist),100)
end