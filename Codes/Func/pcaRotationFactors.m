%
%
% Rotation of factors according to its explained variance
% 
%
% Ziqiang Wei
% version 1.0
%



function newLoadingMat = pcaRotationFactors(loadingMat)
    
    rotationMat    = pcacov(loadingMat' * loadingMat);
    newLoadingMat  = loadingMat * rotationMat;
    
end