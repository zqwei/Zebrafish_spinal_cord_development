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

function CHPoints = smoothedBoundary3(x, y, z, radius)
    pointList     = balledData (x, y, z, radius);
    pindex        = convhull(pointList(:,1), pointList(:,2), pointList(:,3), 'simplify',true);
%     px            = pointList(pindex, 1);
%     py            = pointList(pindex, 2);
%     pz            = pointList(pindex, 3);
%     pxyz          = polygon([px(1:end-1), py(1:end-1), pz(1:end-1)], radius, 30);
    CHPoints      = pointList(pindex, :);
end

function pointList = balledData (x, y, z, radius)
    
    numPoints      = 21 * 21; % matlab default
    pointList      = zeros(length(x)*numPoints, 3);
    [px, py, pz]   = sphere();
    addList        = [px(:), py(:), pz(:)] * radius;
    
    
    for nPoint     = 1:length(x)
        pointList((nPoint-1)*numPoints+(1:numPoints), :) = bsxfun(@plus, [x(nPoint), y(nPoint), z(nPoint)], addList);
    end
    
end


