function tval = tcoef(x,y)
%
% function tval = tcoef(x,y)
% forcing function for finite diff function cdv_2d
%
% tval: rate of reaction/absorption direction
% x   : x-value of point in 2D domain
% y   : y-value of point in 2D domain

% constant absorption
tval = 0.0;

% constant absorption in circle
% if (x - 0.85)^2 + (y - 0.65)^2 <= 0.01
%   tval = 100;
% elseif (x - 0.25)^2 + (y - 0.35)^2 <= 0.01
%   tval = 20;
% else
%   tval = 0;
% end

%

% constant absorption in circle
% if (x - 0.85)^2 + (y - 0.65)^2 <= 0.01
%   tval = 100;
% %elseif (x - 0.25)^2 + (y - 0.35)^2 <= 0.01
% %  tval = 20;
% else
%   tval = 0;
% end

