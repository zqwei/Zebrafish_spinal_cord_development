% 
% Get radius of a cluster of data
%
% Input:
% 
% X        --- N x m
% 
% -------------------------------------------------------------------------
% version 1.0
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 



function r      = getRadius(X)
    
    centriodX   = mean(X, 1);
    distX       = bsxfun(@minus, X, centriodX);
    r           = max(sqrt(sum(abs(distX).^2, 2)));
end