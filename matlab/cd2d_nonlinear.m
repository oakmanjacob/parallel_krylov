function  [F,J,rhsF,Fmat] = cd2d_nonlinear(nx,ny,dx,dy,x0,y0,uvec,...
         p,q,pdx,qdx,r,s,t,f,s_bnd,w_bnd,n_bnd,e_bnd)

% This functions generates a sparse matrix and right hand side from the
% 2-dimensional convection-diffusion equation 
%
%     - (p*u_x)_x - (q*u_y)_y + r*u_x + s*u_y + t*u = f
%
% using 2nd order central finite differences discretization.  Here p and q
% are functions of x,y, and u (and therefore the convection-diffusion
% equation above is nonlinear in diffusion.
%
% nx: number of grid points in x-direction 
%     not counting boundary points for Dir bc
%     counting boundary points for Neumann/Robin bc
% ny: number of grid points in y-direction
%     not counting boundary points for Dir bc
%     counting boundary points for Neumann/Robin bc
% dx: (fixed) mesh width in x-direction; depends on bc
% dy: idem
% x0: x-coeff lower left corner of domain (south west corner)
% y0: y-coeff lower left corner of domain (south west corner)
% u:  solution at the previous Newton iteration, stored as a vector

% The following are all functions of x, y, and u: p(x,y,u), q(x,y,u), etc.
%
% p:  diffusion coefficient (positive) in x-direction
% q:  idem y-direction
% pdx: derivative of p
% qdx: derivative of q
%
% The following are all functions of x and y:
%
% r:  convection in x-direction (flow velocity)
% s:  idem y-direction
% t:  reaction or absorption rate 
% f:  forcing function (rhs); source or sink

% The following are constants (can be modified to be functions of x or y)

% s_bnd: Dirichlet boundary condition - south
% w_bnd: Dirichlet boundary condition - west
% n_bnd: Dirichlet boundary condition - north
% e_bnd: Dirichlet boundary condition - east
%
% Author: Arielle Grim-McNally 2017 (modified from cd_2d.m code, author?)
% Most recent modifications: February 2018

nrows = nx*ny;
dx2 = dx^2;
dy2 = dy^2;

% Allote space for A(u)u and f (to put together -F)
mat_tmpF = zeros(nx,ny,5);
rhs_tmpF = zeros(nx,ny);

% Allote space for J (the Jacobian)
mat_tmpJ = zeros(nx,ny,5);

% Convert the solution vector to a matrix
u_length = length(uvec);
u_size = sqrt(u_length);
u = zeros(u_size);
ucolAdd = zeros(u_size,1);
ucol = 1;
for i = 1:u_length
    idx = mod(i,u_size);
    if idx ~= 0
        ucolAdd(idx) = uvec(i);
    else
        ucolAdd(u_size) = uvec(i);
        u(:,ucol) = ucolAdd;
        ucol = ucol + 1;
        ucolAdd = zeros(u_size,1);
    end
end


for i = 1:nx
  for j = 1:ny
    xc = x0 + i*dx;   yc = y0 + j*dy;   % center
    xs = xc;          ys = yc - 1/2*dy; % south
    xe = xc + 1/2*dx; ye = yc;          % east
    xn = xc;          yn = ys + dy;     % north
    xw = xe - dx;     yw = yc;          % west
    
    % Obtain surrounding u-values: 
    uc = u(i,j); % Center
    
    if j ~= 1 % South
        us = u(i,j-1);
    else
        us = s_bnd(xs,ys);
    end
    
    if i ~= nx % East
        ue = u(i+1,j);
    else
        ue = e_bnd(xe,ye);
    end
    
    if j ~= ny % North
        un = u(i,j+1);
    else
        un = n_bnd(xn,yn);
    end
    
    if i ~= 1 % West
        uw = u(i-1,j);
    else
        uw = w_bnd(xw,yw);
    end
    
    % Approximate the solution at the midpoints by taking the average of
    % the solution at the surrounding points
    south = (uc+us)/2; %q 
    west = (uc+uw)/2; %p
    north = (uc+un)/2; %q
    east = (uc+ue)/2; %p
    
    % diffusion - A(u)
    mat_tmpF(i,j,1) = -q(xs,ys,south)/dy2;  % south
    mat_tmpF(i,j,2) = -p(xw,yw,west)/dx2;  % west
    mat_tmpF(i,j,4) = -p(xe,ye,east)/dx2;  % east
    mat_tmpF(i,j,5) = -q(xn,yn,north)/dy2;  % north
    mat_tmpF(i,j,3) = -(mat_tmpF(i,j,1) + mat_tmpF(i,j,2) + ...
                       mat_tmpF(i,j,4) + mat_tmpF(i,j,5));  % center
    
    % convection - A(u)
    mat_tmpF(i,j,1) = mat_tmpF(i,j,1) - s(xc,yc)/(2*dy); % south
    mat_tmpF(i,j,2) = mat_tmpF(i,j,2) - r(xc,yc)/(2*dx); % west
    mat_tmpF(i,j,4) = mat_tmpF(i,j,4) + r(xc,yc)/(2*dx); % east
    mat_tmpF(i,j,5) = mat_tmpF(i,j,5) + s(xc,yc)/(2*dy); % north
    
    % reaction - A(u)
    mat_tmpF(i,j,3) = mat_tmpF(i,j,3) + t(xc,yc);
    
    % rhs - f
    rhs_tmpF(i,j) = f(xc,yc);
    
    % Build the Jacobian:
    % South (q)
    mat_tmpJ(i,j,1) = 1/dy2*(q(xs,ys,south)+(1/2)*us*qdx(xs,ys,south)-(1/2)*uc*qdx(xs,ys,south));
    
    % West (p)
    mat_tmpJ(i,j,2) = 1/dx2*(p(xw,yw,west)+(1/2)*uw*pdx(xw,yw,west)-(1/2)*uc*pdx(xw,yw,west));
    
    % East (p)
    mat_tmpJ(i,j,4) = 1/dx2*(p(xe,ye,east)-(1/2)*uc*pdx(xe,ye,east)+(1/2)*ue*pdx(xe,ye,east));
    
    % North (q)
    mat_tmpJ(i,j,5) = 1/dy2*(q(xn,yn,north)-(1/2)*uc*qdx(xn,yn,north)+(1/2)*un*qdx(xn,yn,north));
    
    % Center (p and q)
    mat_tmpJ(i,j,3) = 1/(dx*dy)*(-q(xs,ys,south)-p(xw,yw,west)-p(xe,ye,east)-q(xn,yn,north)+...
        (1/2)*us*qdx(xs,ys,south)+(1/2)*uw*pdx(xw,yw,west)+(1/2)*ue*pdx(xe,ye,east)+ (1/2)*un*qdx(xn,yn,north)-...
        uc*((1/2)*qdx(xs,ys,south)+(1/2)*pdx(xw,yw,west)+(1/2)*pdx(xe,ye,east)+(1/2)*qdx(xn,yn,north)));
  end
end

% boundary conditions
xbnd = x0+dx:dx:x0+nx*dx; ybnd = y0+dy:dy:y0+ny*dy;
sbval = s_bnd(xbnd,ybnd)';
ebval = e_bnd(xbnd,ybnd)';
nbval = n_bnd(xbnd,ybnd)';
wbval = w_bnd(xbnd,ybnd)';

% south
rhs_tmpF(1:nx,1) = rhs_tmpF(1:nx,1) - mat_tmpF(1:nx,1,1).*sbval;
mat_tmpF(1:nx,1,1) = 0;
mat_tmpJ(1:nx,1,1) = 0;

% east
rhs_tmpF(nx,1:ny) = rhs_tmpF(nx,1:ny) - mat_tmpF(nx,1:ny,4).*ebval;
mat_tmpF(nx,1:ny,4) = 0;
mat_tmpJ(nx,1:ny,4) = 0;

% north
rhs_tmpF(1:nx,ny) = rhs_tmpF(1:nx,ny) - mat_tmpF(1:nx,ny,5).*nbval;
mat_tmpF(1:nx,ny,5) = 0;
mat_tmpJ(1:nx,ny,5) = 0;

% west
rhs_tmpF(1,1:ny) = rhs_tmpF(1,1:ny) - mat_tmpF(1,1:ny,2).*wbval;
mat_tmpF(1,1:ny,2) = 0;
mat_tmpJ(1,1:ny,2) = 0;

% create sparse matrix and right hand side
rowsF = zeros(5*nrows,1); rowsJ = rowsF;
colsF = zeros(5*nrows,1); colsJ = colsF;
valsF = zeros(5*nrows,1); valsJ = valsF;
rhsF  = zeros(nrows,1);

idx = 0;
for k = 1:nrows
  i = mod(k-1,nx)+1;
  j = (k-i)/nx+1;
  rowsF(idx+1:idx+5) = [k;k;k;k;k];
  rowsJ(idx+1:idx+5) = [k;k;k;k;k];
  
  % For references to boundary points (not in grid/matrix), replace by 
  % incorrect but valid references (in matrix) with zero value. 
  % Those are later automatically removed by "sparse"
  s_idx = max(1,k-nx);
  w_idx = max(1,k-1);
  e_idx = min(k+1,nrows);
  n_idx = min(k+nx,nrows);
  colsF(idx+1:idx+5) = [s_idx;w_idx;k;e_idx;n_idx];
  colsJ = colsF;
  valsF(idx+1:idx+5) = [mat_tmpF(i,j,1); ...
                       mat_tmpF(i,j,2); ...
                       mat_tmpF(i,j,3); ...
                       mat_tmpF(i,j,4); ...
                       mat_tmpF(i,j,5)];
  valsJ(idx+1:idx+5) = [mat_tmpJ(i,j,1); ...
                       mat_tmpJ(i,j,2); ...
                       mat_tmpJ(i,j,3); ...
                       mat_tmpJ(i,j,4); ...
                       mat_tmpJ(i,j,5)];
  rhsF(k) = rhs_tmpF(i,j);
  idx = idx+5;
end

Fmat = sparse(rowsF,colsF,valsF);
J = sparse(rowsJ,colsJ,valsJ);
% Needs to be more efficient....but it works.
% Build the vector F, containing the entries F_(ij) = A(u)u-f
F = Fmat*uvec - rhsF;

end
