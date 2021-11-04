function Mb = rhsLeftPrec(mats,b)
% Function rhsLeftPrec, computes M\b in the case of left preconditioning,
% currently only assumes preconditioner type ILUT(P) & SAMs
% Author:   Arielle (Grim-McNally) Carr, 2017
% Contact:  arg318@lehigh.edu
%
% Inputs:   mats:       Array of matrices, containing ILUTP preconditioner,
%                       factors L and U stored as separate matrices
%           b:          RHS
% Output:   Mb:         M\b

if strcmp('ILU',mats{1}) == 1 %|| strcmp('SAM',mats{1}) == 1% ILU/ILUTP or SAM with ILUTP
    L = mats{3};
    U = mats{4};
    Mb = U\(L\b);
elseif strcmp('SAM',mats{1}) == 1
    L = mats{3};
    U = mats{4};
    M = mats{5};
    Mb = U\(L\(M*b));
elseif strcmp('NONE',mats{1}) == 1
    Mb = b;
end
