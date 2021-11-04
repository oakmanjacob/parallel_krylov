function y = matvec(mats,x,precSide)
% Function matvec.m to determine A*x
% Authors:  Arielle (Grim-McNally) Carr, 2017
% Contact:  arg318@lehigh.edu
%
% Inputs:       mats:       Array of matrices (may include a
%                           preconditioner, at minimum contains the LHS
%                           matrix A), currently only assumes ILUTP and SAM
%                           (with ILUTP) as preconditioner types.
%               x:          vector of length n
%               precSide:   If a preconditioner is specified, this states
%                           which side to apply it from 'Right' or 'Left'
%                           (if not given, default is 'Right')
% Latest update: October 15, 2021
if strcmp('NONE',mats{1}) ==1
    A = mats{2};
    y = A*x;
end
if strcmp('ILU',mats{1}) == 1 % ILU or ILUTP
    A = mats{2};
    L = mats{3};
    U = mats{4};
    if strcmp(precSide,'Right') == 1 || nargin < 3
        y = A*(U\(L\x));
    elseif strcmp(precSide,'Left') == 1
        y = U\(L\(A*x));
    end
end
if strcmp('SAM',mats{1}) == 1 % Sparse Approximate Map with ILUTP
    A = mats{2};
    L = mats{3};
    U = mats{4};
    M = mats{5};
    if strcmp(precSide,'Right') == 1 || nargin < 3
        % y = M*(A*(U\(L\x))); % If mapping from the left
        y = A*(M*(U\(L\x))); % If mapping from the right
    elseif strcmp(precSide,'Left') == 1
        %y = U\(L\(M*(A*x))); % If mapping from the left
        y = U\(L\(A*(M*x))); % If mapping from the right
    end
end


