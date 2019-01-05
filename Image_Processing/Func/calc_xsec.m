%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc_xsec.m
% calculate the midial-lateral (ML) and dorsal-ventral (DV) position of a
% point P based on the cross section plane defined by P, Pa and Pb
% P, Pa, Pb: 3x1 vectors of xyz coordinates
% y: ML position, normalized against range, [-1, 1] midline=0, Pa=-1, Pb=1
% z: DV position, absolute distance of P to PaPb
% width = dist(Pa,Pb)
%
% -------------------------------------------------------------------------
% Yinan Wan
% wany@janelia.hhmi.org
function [y, z, width] = calc_xsec(P, Pa, Pb)
width = norm(Pa - Pb);
z = norm(cross(P-Pa, P-Pb))/norm(Pb-Pa);
y = sqrt(norm(P-Pa)^2 - z^2)*2 - width;
y = y/width;