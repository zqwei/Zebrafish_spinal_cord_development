function [AP, LR, ori] = fitBodyAxes(points, x, y)
% input points: xyz coordinates of points
% x,y,z: AP, LR and DV location in atlas coordinate system
validID = ~isnan(points(:, 1));
points = points(validID, :);
x = x(validID);
y = y(validID);

nPoints = size(points, 1);
ori = nanmean(points);
% linear fit to get AP direction
lx = [ones(nPoints, 1), x]\points;
AP = lx(2, :)/norm(lx(2, :));

% project all points along the AP axis
pointsCrossX = project2Plane(points, AP);
ly = [ones(nPoints, 1), y]\pointsCrossX;
LR = ly(2, :)/norm(ly(2, :));

end
