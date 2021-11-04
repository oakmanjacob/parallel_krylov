function qval = qcoef(x,y,u)
%
% function qval = qcoef(x,y)
% coefficient function for finite diff function cdv_2d
%
% qval: rate of diffusion in y direction
% x   : x-value of point in 2D domain
% y   : y-value of point in 2D domain

% constant diffusion
qval = 0.1+u^2;

% piece-wise constant
% assume [0.1] x [0,1]
% if (y > 0.25) && (y < 0.75) && (x > 0.25) && (x < 0.75),
%   qval = 0.01;
% else
%   qval = 1;
% end
