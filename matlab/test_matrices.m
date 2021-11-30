n = 16384;
x0 = zeros(n, 1);
tol = 10.^-6;

mgmres_times = zeros(1,10);
gmres_times = zeros(1,10);

for k = 1:10
    import = readmatrix(append("../data/n", int2str(n) ,"/A_", int2str(k), ".csv"));
    
    A = sparse(import(:,1),import(:,2),import(:,3), n, n);
    
    b = readmatrix(append("../data/n", int2str(n) ,"/b_", int2str(k), ".csv"));
    
    mats = {"NONE", A};
    t0 = clock;
    [xmGM,rmGM,r_nrmmGM,itermGM,flagmGM] = mgmres(mats,b,x0,tol,n,n, 'Right');
    itermGM
    mgmres_times(k) = round(etime(clock,t0) * 1000);

    t1 = clock;
    [xGM,flagGM,rGM,iterGM,r_nrmGM] = gmres(A,b,n,tol,n);
    iterGM
    gmres_times(k) = round(etime(clock,t1) * 1000);
end

figure(1)
plot(1:10, mgmres_times);
hold on;
plot(1:10, gmres_times);
plot(1:10, [2635]);

legend("mgmres", "gmres", "ogmres")