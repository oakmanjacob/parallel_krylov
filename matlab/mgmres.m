function [x,r,r_nrm,iter,flag] = mgmres(mats,b,x,tol,max_it,m,precSide)
% Function GMRES(m) Restarted GMRES, solves the linear system Ax = b
% Author:           Eric de Sturler <= 2012,
%                   Latest modifications Arielle Carr 2021
% Contact:          arg318@lehigh.edu
%
% Inputs:           mats:       Array of matrices (may include a 
%                               preconditioner, at minimum contains the LHS
%                               matrix A)
%                   b:          RHS 
%                   x:          Initial guess
%                   tol:        Convergence tolerance
%                   max_it:     Maximum number of iterations
%                   m:          Number of iterations before restarting 
%                               GMRES (set to n if using Full GMRES)
%                   precSide:   'Left' or 'Right' preconditioning, if left
%                               blank, will default to Right
% Outputs:          x:      Solution to Ax = b
%                   r:      Residual r = b - Ax
%                   r_nrm:  Residual norm vector
%                   iter:   Number of iterations 
%                   flag:   1: Converged to solution; 
%                           0: Did not converge to solution within max_it

n = length(x);
iter = 0;
r_nrm = zeros(max_it+1,1);
flag = 0;

Ax = matvec(mats,x,precSide);
if strcmp(precSide,'Left') == 1
    Mb = rhsLeftPrec(mats,b);
elseif strcmp(precSide,'Right') == 1 || strcmp('NONE',mats{1}) == 1 || nargin < 7 
    Mb = b;
end
r = Mb - Ax; % Compute the residual
r_nrm(1) = norm( r );
tol = norm(Mb)*tol;

if ( r_nrm(1) <= tol )
    return
end

V(1:n,1:m+1) = zeros(n,m+1);
H(1:m+1,1:m) = zeros(m+1,m);
cs(1:m) = zeros(m,1);
sn(1:m) = zeros(m,1);
e1    = zeros(m+1,1);
e1(1) = 1.0;

iter = 1;
notconv = 1;
while notconv && (iter <= max_it)
    % Arnoldi (Construct orthonormal basis)
    V(:,1) = r / r_nrm(iter);
    s = r_nrm(iter)*e1;       
    i = 1;                    % initialize m steps of Arnoldi
    G = zeros(m+1,m);         % initialize Givens rotation matrix
    while notconv && (iter <= max_it) && (i<=m)
        iter = iter+1;
        w = matvec(mats,V(:,i),precSide);
        for k = 1:i,
            H(k,i)= V(:,k)'*w;
            G(k,i) = H(k,i);
            w = w - H(k,i)*V(:,k);
        end
        H(i+1,i) = norm( w );
        G(i+1,i) = H(i+1,i);
        if H(i+1,i) > 0
            V(:,i+1) = w / H(i+1,i);
            % apply Givens rotation 
            for k = 1:i-1,       
                temp     =  cs(k)*H(k,i) + conj(sn(k))*H(k+1,i);
                H(k+1,i) = -sn(k)*H(k,i) + conj(cs(k))*H(k+1,i);
                H(k,i)   = temp;
            end
            % form i-th Givens rotation matrix
            [cs(i),sn(i)] = rotmat( H(i,i), H(i+1,i) ); 
            
            % approximate residual norm
            temp   = cs(i)*s(i);                        
            s(i+1) = -sn(i)*s(i);
            s(i)   = temp;
            H(i,i) = cs(i)*H(i,i) + conj(sn(i))*H(i+1,i);
            H(i+1,i) = 0.0;
            r_nrm(iter)  = abs(s(i+1));
            
            if ( r_nrm(iter) <= tol ),
                notconv=0;
            else
                i = i+1;
            end
        else
            notconv=0;
        end
    end
    
    % Form the (approximate) solution
    if ( ~notconv ) 
        y = H(1:i,1:i) \ s(1:i);
        x = x + V(:,1:i)*y;
    else
        y = H(1:m,1:m) \ s(1:m);
        x = x + V(:,1:m)*y;
        i = i - 1;
    end
    Ax = matvec(mats,x,precSide);
    if strcmp(precSide,'Left') == 1
         Mb = rhsLeftPrec(mats,b);
    elseif strcmp(precSide,'Right') == 1 || strcmp('NONE',mats{1}) == 1 || nargin < 7
        Mb = b;
    end
    r = Mb - Ax; % Compute residual
    
    r_nrm(iter) = norm(r); % Compute norm of residual
    
    if (r_nrm(iter) <= tol), % check convergence
        notconv = 0;
    else
        notconv = 1;
    end
end
if ( r_nrm(iter) <= tol ) % converged
    flag = 1; 
end
% eliminate superfluous zeros in vector of residuals norms
r_nrm = r_nrm(1:iter); 
end

