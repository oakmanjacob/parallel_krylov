function fval = fcoef(x,y)
%
% function fval = fcoef(x,y)
% forcing function for finite diff function cdv_2d
%
% fval: rate of convection in y direction
% x   : x-value of point in 2D domain
% y   : y-value of point in 2D domain

% constant forcing function
% fval = 0;

% small sources on line
fval = 0;%cos(x)+sin(y);
% if (0.005<y) && (y<0.015),
%   if (0.415<x) && (x<=0.425),
%     fval = 1;
%   elseif (0.29<x) && (x<0.31),
%     fval = 0;
%   elseif (0.39<x) && (x<0.41),
%     fval = 1;
%   elseif (0.49<x) && (x<0.51),
%     fval = 1;
%   elseif (0.59<x) && (x<0.61),
%     fval = 0;
%   elseif (0.69<x) && (x<0.71),
%     fval = 1;
%   else
%     fval = 0;
%   end
% end
      