function [x, y, z, width, a, b] = convert2atlas3D(points, ra, rb)
% convert the points into 3D atlas based on reference curves a and b
% points: nPoints x 3 matrix, list of points to be converted
% a, b: N x 3 matrix containing the reference points in order
% x: m x 1 vector of the AP position by segment
% y, z, width as in calc_xsec
sampling = 100;
[ida, a] = interp_curve(ra, sampling);
[idb, b] = interp_curve(rb, sampling);


nPoints = size(points, 1);
x = zeros(nPoints, 1);
y = x; z = x; width = x;
for i = 1:nPoints
    P = points(i, :);
    % for each point to be processed, find the closest point on each curve
    [da, iPa] = min(sqrt(sum(abs(repmat(P, numel(ida), 1)- a).^2,2)));    
    [db, iPb] = min(sqrt(sum(abs(repmat(P, numel(idb), 1)- b).^2,2)));
    
    if da <= db
            x(i) = ida(iPa);
            Pa = a(iPa, :);
            Pb = b(iPa, :);
    else
            x(i) = ida(iPb);
            Pb = b(iPb, :);
            Pa = a(iPb, :);
    
    end
    [y(i), z(i), width(i)] = calc_xsec(P, Pa, Pb);
end
    
    