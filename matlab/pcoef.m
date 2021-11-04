function pval = pcoef(x,y,u)
%
% function pval = pcoef(x,y)
% coefficient function for finite diff function cdv_2d
%
% pval: rate of diffusion in x direction
% x   : x-value of point in 2D domain
% y   : y-value of point in 2D domain

% constant diffusion
pval = 0.1+u^2;

% piece-wise constant
% assume [0,1] x [0,1]
% if (y > 0.25) && (y < 0.75) && (x > 0.25) && (x < 0.75),
%   pval = 0.01;
% else
%   pval = 1;
% end
