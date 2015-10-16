% 
% Plot circle at location x, y, with radius r and color c
%
% Input:
% 
% X        --- N x 3
% 
% -------------------------------------------------------------------------
% version 1.0
%
% 
% -------------------------------------------------------------------------
% Ziqiang Wei
% weiz@janelia.hhmi.org
% 



function drawCircle(x, y, r, c)
    t           = linspace(0,2*pi);
    hold on
    for nCircle = 1:length(x)
        if size(c,1)== 1
            fill(x(nCircle)+r(nCircle)*cos(t), y(nCircle)+r(nCircle)*sin(t), c, 'edgecolor', c);
        else
            fill(x(nCircle)+r(nCircle)*cos(t), y(nCircle)+r(nCircle)*sin(t), c(nCircle,:), 'edgecolor', c(nCircle,:));
        end
    end
    hold off
end