%% Fichier de test comparant les differents solveurs selon la taille de la fenetre


clear all;
close all;
load mat1.mat;

tol = 1e-10;
n = size(A, 1);
maxit = 2*n;
x0 = zeros(n,1);

[M1,M2] = ilu(A);

figure(1)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,1,1)


% test avec m = 10
m = 10;
normRHS = norm(M2\(M1\b));

[x, flag, resrel, iter, resvec] = restartedGMRES(A, b, m, tol, maxit, M1, M2, x0);
semilogy(resvec/normRHS);
hold on



[x, flag, resrel, iter, resvec] = dqGMRES1(A, b, m, tol, maxit, M1, M2, x0);
semilogy(resvec/normRHS,'o');
hold on

[x, flag, resrel, iter, resvec] = dqGMRES2(A, b, m, tol, maxit, M1, M2, x0);
semilogy(resvec/normRHS);
hold on


title(['mat1, m = 10'])
legend('restartedGMRES','dqGMRES1','dqGMRES2');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,1,2)


% test avec m = 20
m = 20;
normRHS = norm(M2\(M1\b));

[x, flag, resrel, iter, resvec] = restartedGMRES(A, b, m, tol, maxit, M1, M2, x0);
semilogy(resvec/normRHS);
hold on



[x, flag, resrel, iter, resvec] = dqGMRES1(A, b, m, tol, maxit, M1, M2, x0);
semilogy(resvec/normRHS,'o');
hold on

[x, flag, resrel, iter, resvec] = dqGMRES2(A, b, m, tol, maxit, M1, M2, x0);
semilogy(resvec/normRHS);
hold on


title(['mat1, m = 20'])
legend('restartedGMRES','dqGMRES1','dqGMRES2');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,1,1)


% test avec m = 30
m = 30;
normRHS = norm(M2\(M1\b));

[x, flag, resrel, iter, resvec] = restartedGMRES(A, b, m, tol, maxit, M1, M2, x0);
semilogy(resvec/normRHS);
hold on



[x, flag, resrel, iter, resvec] = dqGMRES1(A, b, m, tol, maxit, M1, M2, x0);
semilogy(resvec/normRHS,'o');
hold on

[x, flag, resrel, iter, resvec] = dqGMRES2(A, b, m, tol, maxit, M1, M2, x0);
semilogy(resvec/normRHS);
hold on


title(['mat1, m = 30'])
legend('restartedGMRES','dqGMRES1','dqGMRES2');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,1,2)


% test avec m = 50
m = 50;
normRHS = norm(M2\(M1\b));

[x, flag, resrel, iter, resvec] = restartedGMRES(A, b, m, tol, maxit, M1, M2, x0);
semilogy(resvec/normRHS);
hold on


[x, flag, resrel, iter, resvec] = dqGMRES1(A, b, m, tol, maxit, M1, M2, x0);
semilogy(resvec/normRHS,'o');
hold on

[x, flag, resrel, iter, resvec] = dqGMRES2(A, b, m, tol, maxit, M1, M2, x0);
semilogy(resvec/normRHS);
hold on


title(['mat1, m = 50'])
legend('restartedGMRES','dqGMRES1','dqGMRES2');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%