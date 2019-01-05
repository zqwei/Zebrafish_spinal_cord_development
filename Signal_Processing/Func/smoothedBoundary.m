%
%
% smoothed convex hull
% 
%
% Ziqiang Wei
% version 1.0
%
%

function CHPoints = smoothedBoundary(x, y, radius)
    pointList     = balledData (x, y, radius/2);
    shrinkRatio   = 0.1; % 0 for top view
    pindex        = boundary(pointList(:,1), pointList(:,2), shrinkRatio);
    px            = pointList(pindex, 1);
    py            = pointList(pindex, 2);
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


