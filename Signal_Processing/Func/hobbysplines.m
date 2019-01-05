%HOBBYSPLINES Draws open or closed smooth curves from keypoints
%
% hobbysplines({[x1 y1], [x2 y2],... ,[xN yN]},[opts])
% hobbysplines({[x1 y1 z1], [x2 y2 z2],... ,[xN yN zN]},[opts])
%
% Draws a closed (cyclic) curve through keypoints 1 to N.
% Keypoints may be specified with optional slopes and tension parameters.
% The syntax to do this (replacing `[x1 y1]` above, say) is
%
% 2D:
%    { [x1 y1] s1 t1 }
%    { [x1 y1] s1 t1_in t1_out }
%    { [x1 y1] [u1 v1] ... }
%
% 3D:
%    { [x1 y1 z1] [u1 v1 w1] t1 }
%    { [x1 y1 z1] [u1 v1 w1] t1_in t1_out }
%
% where `s1` is the slope in degrees of the curve through `[x1 y1]`,
% `[u1 v1 w1]` is the unit vector describing the slope of the curve,
% and `t1_out` & `t1_in` are the "exit" and "entry" tensions of the curve
% approaching that point. If `t1` only is specified, this is both the
% entry and the exit tension.
%
% Use '' to indicate default values here; for example, for a default slope
% but specified "entry" tension of 2.0, use
%
%    { [x1 y1] '' 2.0 ''}
%
% and note that trailing optional arguments can also be omitted.
% According to Hobby, tensions should not be specified less than 0.75.
% This software does not enforce this restriction.
% Note that a tension of 1 creates approximately circular plots.
%
% Optional arguments given by [opts] can be any combination of the
% following:
%
%   OPTION       DEFAULT            DESCRIPTION
%   ------       -------            -----------
%   'tension'    [1]                default tension between points
%   'offset'     [0 0 0]            offset to add to each control point
%   'cycle'      [true]             whether to draw a cyclic curve
%   'debug'      [false]            draw and label keypoints on the curve
%   'linestyle'  {'linewidth',1}    line style option(s)
%   'color'      'black'            colour of the curve
%
% Distributed under the terms and conditions of the 2-clause BSD license:  
% <http://opensource.org/licenses/bsd-license.php>
%
% Copyright 2013 Will Robertson and The University of Adelaide
% All rights reserved.
%
%
% -----------------------------------------------------------------------
% Input:
% x, y        ---- location of polygon
% px, py      ---- smoothed points
%
%
%
% Modified by Ziqiang Wei
% weiz@janelia.hhmi.org
%



function pxyz         = hobbysplines(points)
    
    [numPoints, numDim] = size(points);
    
    
    % computing slope of curve ----
    if numDim         < 3
        points            = [points, zeros(numPoints, 1)];
    end
    points(end+1,:)   = points(1,:);
    pointDirs         = zeros(numPoints+1, 3);
    
    pxyz              = [];
    
    for nDir          = 1:numPoints+1
        switch nDir
            case 1
                pointDirs(nDir, :) = points(2, :) - points(end-1, :);
            case numPoints+1
                pointDirs(nDir, :) = points(2, :) - points(end-1, :);
            otherwise
                pointDirs(nDir, :) = points(nDir+1, :) - points(nDir-1, :);
        end
    end
    
    for nPoint        = 1:numPoints
        theta         = arg(pointDirs(nPoint, :))-arg(points(nPoint+1, :) - points(nPoint, :));
        phi           = arg(points(nPoint+1, :) - points(nPoint, :))-arg(pointDirs(nPoint+1, :));
        [rho,sigma]   = velocity_parameters(theta,phi);
        P1            = points(nPoint, :);
        P2            = points(nPoint, :)  +rho/3  *norm(points(nPoint+1, :) - points(nPoint, :))*pointDirs(nPoint, :);
        P3            = points(nPoint+1, :)-sigma/3*norm(points(nPoint+1, :) - points(nPoint, :))*pointDirs(nPoint+1, :);
        P4            = points(nPoint+1, :);
        txyz          = plot_bezier(P1,P2,P3,P4);
        pxyz          = [pxyz; txyz]; %#ok<*AGROW>
    end
end

function o = arg(w)
  o = atan2(w(2),w(1));
end

function [rho,sigma] = velocity_parameters(theta,phi)
    % From "Smooth, easy to compute interpolating splines" by John D. Hobby
    % <http://www.springerlink.com/content/p4w1k8w738035w80/>
    a = 1.597;
    b = 0.07;
    c = 0.37;

    st = sin(theta);
    ct = cos(theta);
    sp = sin(phi);
    cp = cos(phi);

    alpha = a*(st-b*sp)*(sp-b*st)*(ct-cp);
    rho   = (2+alpha)/(1+(1-c)*ct+c*cp);
    sigma = (2-alpha)/(1+(1-c)*cp+c*ct);

end

function Q = plot_bezier(P1,P2,P3,P4)
    N = 50;
    t = linspace(0,1,N)';
    c1 = 3*(P2 - P1);
    c2 = 3*(P1 - 2*P2 + P3);
    c3 = -P1 + 3*(P2-P3) + P4;
    Q = t.^3*c3 + t.^2*c2 + t*c1 + repmat(P1,[N 1]);
end
