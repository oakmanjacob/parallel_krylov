function [p,back] = backstep(back,u,nx,ny,dx,dy,x0,y0,max_m,F,p,alpha)
normF = norm(F);
u_prop = u + p;
[F_prop,~,~,~] = cd2d_nonlinear(nx,ny,dx,dy,x0,y0,u_prop, ...
    @pcoef,@qcoef,@pcoefdx,@qcoefdx,@rcoef,@scoef,@tcoef,@fcoef, ...
    @sbc,@wbc,@nbc,@ebc);

m = 0;
half = m;
while (m <=max_m) && (norm(F_prop) > (1-alpha*half)*normF)
    m = m + 1;
    half = 2^(-m);
    p = half*p;
    u_prop = u + p;
    [F_prop,~,~,~] = cd2d_nonlinear(nx,ny,dx,dy,x0,y0,u_prop, ...
        @pcoef,@qcoef,@pcoefdx,@qcoefdx,@rcoef,@scoef,@tcoef,@fcoef, ...
        @sbc,@wbc,@nbc,@ebc);
    back = back + 1;
end

end