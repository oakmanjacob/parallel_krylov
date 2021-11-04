function [ c, s ] = rotmat( a, b )

%
% Compute the Givens rotation matrix parameters for a and b.
%
   if ( b == 0.0 ),
      c = 1.0;
      s = 0.0;
   elseif ( abs(b) > abs(a) ),
      temp = abs(a) / abs(b);
      c = sqrt( temp^2/(temp^2 + 1) );
      s = conj(a/(temp*b)) / sqrt(1 + temp^2);
   else
      temp = abs(b) / abs(a);
      c = 1.0 / sqrt(1 + temp^2);
      s = sqrt(temp^2/(temp^2+1)) * (b/(a*temp));
   end
