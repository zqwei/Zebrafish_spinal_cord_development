% project a point cloud v (N x 3) onto a plane with normal vector n (1 x 3)
function w = project2Plane(v, n)
u = v;
for i = 1:size(v, 1)
    u(i, :) = dot(v(i, :), n)/(norm(n)^2)* n;
end
w = v - u;
end