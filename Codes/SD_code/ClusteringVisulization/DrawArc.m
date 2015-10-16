function DrawArc(arcStartPoint,arcEndPoint,arcCenterPoint)
 v1 = arcStartPoint-arcCenterPoint;
 v2 = arcEndPoint-arcCenterPoint;
 v3 = [0 -1;1 0]*v1;
 a = linspace(0,mod(atan2(det([v1,v2]),dot(v1,v2)),2*pi));
 v = v1*cos(a)+v3*sin(a);
 plot(v(1,:)+arcCenterPoint(1),v(2,:)+arcCenterPoint(2),'k')
 axis equal
end