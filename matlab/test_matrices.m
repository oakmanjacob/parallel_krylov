clear all;
clc;
n = 4096;
x0 = zeros(n, 1);
tol = 10.^-6;

mgmres_times = zeros(1,10);
gmres_times = zeros(1,10);

restart = 1000;

for k = 4:4
    import = readmatrix(append("../data/n", int2str(n) ,"/A_", int2str(k), ".csv"));
    
    A = sparse(import(:,1),import(:,2),import(:,3), n, n);
    
    b = readmatrix(append("../data/n", int2str(n) ,"/b_", int2str(k), ".csv"));
    
    mats = {"NONE", A};
    t0 = clock;
    [xmGM,rmGM,r_nrmmGM,itermGM,flagmGM] = mgmres(mats,b,x0,tol,n,restart, 'Right');
    mgmres_times(k) = round(etime(clock,t0) * 1000);
    flagmGM
    itermGM
    norm(rmGM)
    
    clear A;
    clear b;
    clear import;
    
    import = readmatrix(append("../data/n", int2str(n) ,"/A_", int2str(k), ".csv"));
    
    A = sparse(import(:,1),import(:,2),import(:,3), n, n);
    
    b = readmatrix(append("../data/n", int2str(n) ,"/b_", int2str(k), ".csv"));

    t1 = clock;
    [xGM,flagGM,rGM,iterGM,r_nrmGM] = gmres(A,b,restart,tol,n);
    gmres_times(k) = round(etime(clock,t1) * 1000);
    flagGM
    iterGM
    norm(rGM)
    clear A;
    clear b;
    clear import;
end

figure(1)
plot(1:10, mgmres_times);
hold on;
plot(1:10, gmres_times);

legend("mgmres", "gmres")
hold off;