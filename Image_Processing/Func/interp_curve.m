%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interp_curve.m
% interpolate a 3D curve by sampling times, and return the index range
% curve: N x 3: xyz coordinates of the control points on the curve
%               default index: 1 to N
% index: linspace(-1, N+2, sampling*N) index of each points on the interp curve
% new_curve: sampling*N x 3 matrix: xyz coordinates of new curve points
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
function [index, new_curve] = interp_curve(curve, sampling)

N = size(curve, 1);
index = linspace(0, N+1, N*sampling);

% using cubic spline
% pp = spline(1:N, curve');
% new_curve = ppval(pp, index);
% new_curve = new_curve';

% using 1d pchip interpolation with extrapolation option
% new_curve = interp1(1:N, curve, index, 'pchip', 'extrap');

% quadratic function fit
px = polyfit(1:N, curve(:, 1)', 2);
py = polyfit(1:N, curve(:, 2)', 2);
pz = polyfit(1:N, curve(:, 3)', 2);
new_curve = [polyval(px, index); polyval(py, index); polyval(pz, index)];
new_curve = new_curve';