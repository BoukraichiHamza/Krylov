%% Fichier de test pour valider les différents solveur

clear all;
close all;
load mat1.mat;

tol = 1e-10;
n = size(A, 1);
maxit = n;
x0 = zeros(n,1);
m = 30;

[M1,M2] = ilu(A);

normRHS = norm(M2\(M1\b));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
[x, flag, resrel, iter, resvec] = gmres(A, b, [], tol, maxit, M1, M2, x0);

semilogy(resvec/normRHS, 'x');
hold on;

m = 30;
[x, flag, resrel, iter, resvec] = pGMRES(A, b, 0, tol, maxit, M1, M2, x0);
semilogy(resvec/normRHS, 'o');
title('Validation pGMRES');
legend('gmres Matlab','pGMRES');
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)
[x, flag, resrel, iter, resvec] = gmres(A, b, m, tol, maxit/m, M1, M2, x0);

semilogy(resvec/normRHS, 'x');
hold on;

[x, flag, resrel, iter, resvec] = restartedGMRES(A, b, m, tol, maxit, M1, M2, x0);
semilogy(resvec/normRHS, 'o');
title('Validation restartedGMRES');
legend('gmres Matlab','restartedGMRES');
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3)
[x, flag, resrel, iter, resvec] = gmres(A, b, m, tol, maxit, M1, M2, x0);

semilogy(resvec/normRHS, 'x');
hold on;


[x, flag, resrel, iter, resvec] = dqGMRES1(A, b, m, tol, maxit, M1, M2, x0);
semilogy(resvec/normRHS, 'o');
title('Validation dqGMRES1');
legend('gmres Matlab','dqGMRES1');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4)
[x, flag, resrel, iter, resvec] = gmres(A, b, m, tol, maxit, M1, M2, x0);

semilogy(resvec/normRHS, 'x');
hold on;

[x, flag, resrel, iter, resvec] = dqGMRES2(A, b, m, tol, maxit, M1, M2, x0);
semilogy(resvec/normRHS, 'o');
title('Validation dqGMRES2');
legend('gmres Matlab','dqGMRES2');