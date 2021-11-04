function sval = sbc(x,y)
%
% function sval = sbc(x,y)
%
% implements south Dirichlet boundary condition

% constant

sval = 0.2+y.*(1-y.^2);