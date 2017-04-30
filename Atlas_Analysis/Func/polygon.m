function  xy = polygon(xy,r,res);

% POLYGON  connect Points to Polygon with round edges
%
% XYP = POLYGON( XY , Radius , Resolution );
%
% XY contains the X and Y - Coordinates of the N Points
%  in 2 Columns!
% The first and last Point will be connected (N Edges).
%
% XYP contains the Coordinates of the round edges in 2 Columns
%  size(XYP,1) == Resolution*N
%
% default: Resolution = 10;
%
% Radius could be a Vector of the same Size of XY, then
%  the first Value is valid for the Edge at Point 1, the second at Point 2, ...
%  the last at the last Point.
%
% example:
%
% % 1. Triangle
% xy = [ 2  0 
%        1  1.5
%        0   0   ];
%
% xy = polygon(xy,0.25);
%
% fill(xy(:,1),xy(:,2),[1 0 0])
%
% set(gca,'dataaspectratio',[1 1 1])
%
% % 2. Rectangle
% xy = [ 0.1  0.1 
%        0.9  0.1
%        0.9  0.7
%        0.1  0.7   ];
%
% xy = polygon(xy,[0.1 0.2 0.2 0.0]);
%
% fill(xy(:,1),xy(:,2),[1 0 0])
%
% set(gca,'dataaspectratio',[1 1 1], ...
%         'xlim',[0 1],'ylim',[0 1])
%
% 

if nargin < 3
 res = 10;
end

if ( size(xy,1) == 2  )  &  ( size(xy,2) ~= 2 ) 
 xy = xy';
end

N = size(xy,1);

xy = xy([ (1:N) 1  2 ],:);

ind = ( 1 : size(xy,1)-1 );

v = xy(ind+1,:)-xy(ind+0,:);
  
b = sum(v.*v,2) .^ (0.5);


N   =  size(v,1) - 1;
ind = ( 1 : N );

% Normale on first Line of Edge
n = v(ind,[2 1]) .* ( ones(N,1) * [ 1 -1 ] );
n = n ./ ( b(ind) * [ 1  1 ] );

s = sum(v(ind,:).*v(ind+1,:),2);                      % Scalar
k = v(ind,1).*v(ind+1,2)-v(ind,2).*v(ind+1,1);  % CrossZ

% k  > 0   Rotation positiv
%    < 0            negativ

% Show Inside the Angle  
n = n .* (-sign(k) * [ 1  1 ] );

a = acos(abs(s)./b(ind)./b(ind+1));

% a = asin(abs(k)./b(ind)./b(ind+1))*180/pi;

a = a + ( -a + (pi-a) ) .* ( sign(s) > 0 );


  r = r(:);
if size(r,1) == 1
  r = r*ones(N,1);
elseif size(r,1) < N 
  r = [ r ; zeros(N-size(r,1),1) ];
end
  r = r([ (2:N) 1 ]);


% TangentSegment

t = abs( r ./ sin(a/2) .* sin(pi/2-a/2) );

if any( t+t([(2:N) 1]) > b(ind)  )
%   fprintf([ 'Warning: Radius too large.' char(10)] )
end

% TangentPoint of first Line of Edge
p = xy(ind+1,:) - ...
     ( t * [ 1  1 ] ) .* v(ind,:) ./ ( b(ind) * [ 1  1 ] );

% p2 = xy(ind+1,:) + ...
%      ( t * [ 1  1 ] ) .* v(ind+1,:) ./ ( b(ind+1) * [ 1  1 ] );

% Centre of Radius
m = p + ( r * [ 1  1 ] ) .* n;

% StartAngle of each Circle
a0 = atan2(p(:,2)-m(:,2),p(:,1)-m(:,1));

% AngleRotation of Circle
da = sign(k) .* (pi-a);


ind                    = zeros(res*N,1);
ind(1:res:res*(N-1)+1) = 1;
ind                    = cumsum(ind,1);


% Factor for AngleIncrement

f          = zeros(res*N,1);
f(2:res*N) = 1/(res-1);
f          = cumsum(f);
f          = f-f(ind*res)+1;

ang = a0(ind) + f.*da(ind);

xy = m(ind,:) + ( r(ind)*[ 1  1 ] ) .* [ cos(ang) sin(ang) ];  
