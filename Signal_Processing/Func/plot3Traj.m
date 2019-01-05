% 
% Plot 3D trajectory
%
% Input:
% 
% X        --- N x 3
% 
% -------------------------------------------------------------------------
% version 1.0
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 



function plot3Traj(X, varargin)
    if size(X,2) ~=3
        error('Input needs to be N x 3 matrix; N is the number of observations.')
    else
        plot3(X(:,1),X(:,2),X(:,3), varargin{:});
    end
end