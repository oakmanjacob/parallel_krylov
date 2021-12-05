clear all;
clc;
n = 16384;
x0 = zeros(n, 1);
tol = 10.^-6;

mgmres_times = zeros(1,10);
gmres_times = zeros(1,10);

restart = 1000;

for k = 1:10
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

ogmres = [1544, 6720, 5598, 4314, 5302, 6296, 6116, 5680, 5796, 6983];
ogmres_simd = [1031, 4147, 3940, 3598, 3688, 3965, 4122, 3820, 3697, 3914];

bar([mgmres_times; gmres_times; ogmres; ogmres_simd].');
hold on;

legend("mgmres", "gmres", "ogmres", "ogmres simd")
ylabel("Time (ms)")
xlabel("Convection-Diffusion System")

title("Execution Time of " + n + "x" + n + " C-D Systems")
hold off;

%%

for k = 1:10
    figure(k)

    y = [[1279, 4286, 3009, 3325, 3556, 5116, 4538, 3356, 4093, 3891];[813, 3457, 2137, 1419, 2691, 3375, 2914, 2261, 2242, 2630];[679, 2605, 1423, 843, 1832, 2450, 2238, 1846, 1577, 2242];[641, 2179, 1311, 1306, 1582, 2350, 2049, 1649, 1885, 1821];[644, 2442, 1087, 1361, 1535, 2392, 1832, 1619, 1773, 1766];[624, 2141, 1328, 1410, 1746, 1643, 1899, 1862, 1551, 1725];[647, 2226, 1373, 1313, 1361, 1897, 1949, 1472, 1740, 1410];[654, 2010, 1318, 1318, 1448, 2000, 1840, 1461, 1525, 1742];[1197, 2146, 1556, 1483, 1558, 2043, 1663, 1525, 1462, 1671];[655, 2039, 1389, 1680, 1618, 2497, 1925, 1753, 1644, 2216]];
    plot([10,20,30,40,50,60,70,80,90,100], y(:,k))
    hold on;

    x = [[586, 7631, 2634, 4321, 5243, 7488, 4621, 4555, 6112, 4325];[684, 1110, 1177, 1251, 774, 1729, 1322, 1458, 1412, 1442];[493, 1893, 1026, 761, 787, 1242, 964, 960, 1017, 794];[373, 1026, 842, 1036, 833, 995, 740, 956, 938, 662];[319, 1278, 650, 1069, 758, 1124, 947, 822, 950, 796];[342, 1121, 837, 816, 803, 1048, 930, 993, 951, 911];[610, 1127, 680, 851, 739, 1088, 1164, 771, 1049, 998];[599, 1263, 758, 784, 824, 1129, 948, 827, 1351, 1015];[589, 1380, 817, 995, 888, 1178, 1067, 959, 1325, 1465];[759, 1472, 965, 843, 1022, 1519, 1101, 1191, 1369, 1425];[759, 1958, 1350, 884, 1720, 2290, 1357, 1641, 1567, 1584]];
    plot([4,10,20,30,40,50,60,70,80,90,100], x(:,k))
    hold off;
end