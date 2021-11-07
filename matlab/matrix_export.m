nx = 10;
ny = 10;
dx = 1/(nx+1);
dy = 1/(ny+1);
x0 = 0;
y0 = 0;
n = nx*ny;
u = ones(n,1);

[F,Jac,~,~] = cd2d_nonlinear(nx,ny,dx,dy,x0,y0,u, ...
        @pcoef,@qcoef,@pcoefdx,@qcoefdx,@rcoef,@scoef,@tcoef,@fcoef, ...
        @sbc,@wbc,@nbc,@ebc);
    
[row col v] = find(Jac)
writematrix([row col v], "../data/Jac.csv")
writematrix(F,"../data/F.csv")
    