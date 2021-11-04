
%Computes an ILUTP preconditioner for each system in the sequence
function [iterations, soltime, prectime,fdtime,gmresinfo,back,sol] = ILUTPall_gmres(nx,ny,n,u0)
dx = 1/(nx+1);
dy = 1/(ny+1);
x0 = 0;
y0 = 0;


stopIter = 1000;
iterations = zeros(stopIter,1);
soltime = zeros(stopIter,1);
prectime = zeros(stopIter,1);
fdtime = cell(stopIter,2);
back = 0;
gmresinfo = cell(stopIter,2);
sol = cell(stopIter,1);
u = u0;
ftol = 1.e-4;
F = inf;
p = inf;
k = 0;
alpha = 1.e-4;
max_m = 10;
rtol = 5.e-10;
atol = 5.e-10;
droptol = 1.e-5; lfil = 10;
xx0 = zeros(n,1); ttol = 1.e-10;  max_it = n; restart = 50;

while  (norm(F) > (ftol*rtol + atol)) && (k <= stopIter) %&&norm(p) > ptol %
    k = k + 1
    norm(F)
    tic
    [F,Jac,~,~] = cd2d_nonlinear(nx,ny,dx,dy,x0,y0,u, ...
        @pcoef,@qcoef,@pcoefdx,@qcoefdx,@rcoef,@scoef,@tcoef,@fcoef, ...
        @sbc,@wbc,@nbc,@ebc);
    fdtime{k,1} = toc;
    if k == 1
        ftol = norm(F);
    end
    tic
    [L,U,~,~,Jp] = my_ILUTP(Jac,droptol,lfil);
    prectime(k) = toc;
    mats = {'ILU',Jp,L,U};

    tic
    [vtemp,r,r_nrm,iter,~] = mgmres(mats,F,xx0,ttol,max_it,restart,'Right');
    p = (U\(L\vtemp)); 
    soltime(k) = toc;
    iterations(k) = iter;
    gmresinfo{k,1} = r;
    gmresinfo{k,2} = r_nrm;  
    tic
    [p,back] = backstep(back,u,nx,ny,dx,dy,x0,y0,max_m,F,p,alpha);
    fdtime{k,2} = toc;
    u = u + p;   
    sol{k} = u;
end
