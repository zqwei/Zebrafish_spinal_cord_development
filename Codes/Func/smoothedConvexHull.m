%
%
% smoothed convex hull
% 
%
% Ziqiang Wei
% version 1.0
%
%
%

function CHPoints = smoothedConvexHull(x, y, radius)
    pointList     = balledData (x, y, radius);
    pindex        = convhull(pointList(:,1), pointList(:,2));
    px            = pointList(pindex, 1);
    py            = pointList(pindex, 2);
%     pxyz          = hobbysplines([px, py]);
    pxyz          = polygon([px(1:end-1), py(1:end-1)], radius, 30);
    CHPoints      = pxyz(:, 1:2);
end

function pointList = balledData (x, y, radius)
    
    numPoints      = 36;
    addPoints      = linspace(-pi, pi, numPoints);
    addList        = [cos(addPoints'), sin(addPoints')] * radius;
    
    pointList      = zeros(length(x)*36, 2);
    
    
    for nPoint     = 1:length(x)
        pointList((nPoint-1)*numPoints+(1:numPoints), :) = bsxfun(@plus, [x(nPoint) y(nPoint)], addList);
    end
    
end
