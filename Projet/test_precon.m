%% Fichier de test comparant les differents préconditionneurs sur le solveur pGMRES

clear all;
close all;
load mat1.mat;

tol = 1e-10;
n = size(A, 1);
maxit = 2*n;
x0 = zeros(n,1);

[M11,M21] = ilu(A);

M12 = diag(diag(A));
M22 = eye(n);

M13 = eye(n);
M23 = M13';

M14 = ichol(A);
M24 = M14';
figure(1)
subplot(2,1,1)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test avec m > iter == gmres

normRHS = norm(M21\(M11\b));
m = 30;
[x, flag, resrel, iter, resvec] = pGMRES(A, b, 0, tol, maxit, M11, M21, x0);
semilogy(resvec/normRHS,'o');
hold on

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

normRHS = norm(M22\(M12\b));
m = 30;
[x, flag, resrel, iter, resvec] = pGMRES(A, b, 0, tol, maxit, M12, M22, x0);
semilogy(resvec/normRHS);
hold on
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

normRHS = norm(M23\(M13\b));
m = 30;
[x, flag, resrel, iter, resvec] = pGMRES(A, b, 0, tol, maxit, M13, M23, x0);
semilogy(resvec/normRHS);
hold on

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

normRHS = norm(M24\(M14\b));
m = 30;
[x, flag, resrel, iter, resvec] = pGMRES(A, b, 0, tol, maxit, M14, M24, x0);
semilogy(resvec/normRHS);
title(['mat1, preconditionneur différent'])
legend('ilu','diag','none','ichol');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ A, matname] = bfw398a;
n = size(A, 1);
b = (1:n)';
maxit = 2*n;
x0 = zeros(n,1);

[M11,M21] = ilu(A);

M12 = diag(diag(A));
M22 = eye(n);

M13 = eye(n);
M23 = M13';

M14 = ichol(A);
M24 = M14';
subplot(2,1,2)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test avec m > iter == gmres

normRHS = norm(M21\(M11\b));
m = 30;
[x, flag, resrel, iter, resvec] = pGMRES(A, b, 0, tol, maxit, M11, M21, x0);
semilogy(resvec/normRHS);
hold on

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

normRHS = norm(M22\(M12\b));
m = 30;
[x, flag, resrel, iter, resvec] = pGMRES(A, b, 0, tol, maxit, M12, M22, x0);
semilogy(resvec/normRHS);
hold on
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

normRHS = norm(M23\(M13\b));
m = 30;
[x, flag, resrel, iter, resvec] = pGMRES(A, b, 0, tol, maxit, M13, M23, x0);
semilogy(resvec/normRHS);
hold on

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

normRHS = norm(M24\(M14\b));
m = 30;
[x, flag, resrel, iter, resvec] = pGMRES(A, b, 0, tol, maxit, M14, M24, x0);
semilogy(resvec/normRHS);
title(['BFW398A, preconditionneur différent'])
legend('ilu','diag','none','ichol');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ A, matname] = gre_1107;
n = size(A, 1);
b = (1:n)';
maxit = 2*n;
x0 = zeros(n,1);

[M11,M21] = ilu(A);

M12 = diag(diag(A));
M22 = eye(n);

M13 = eye(n);
M23 = M13';

figure(2)
subplot(2,1,1)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test avec m > iter == gmres

normRHS = norm(M21\(M11\b));
m = 30;
[x, flag, resrel, iter, resvec] = pGMRES(A, b, 0, tol, maxit, M11, M21, x0);
semilogy(resvec/normRHS,'x');
hold on

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

normRHS = norm(M22\(M12\b));
m = 30;
[x, flag, resrel, iter, resvec] = pGMRES(A, b, 0, tol, maxit, M12, M22, x0);
semilogy(resvec/normRHS);
hold on
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

normRHS = norm(M23\(M13\b));
m = 30;
[x, flag, resrel, iter, resvec] = pGMRES(A, b, 0, tol, maxit, M13, M23, x0);
semilogy(resvec/normRHS);
hold on

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

title(['GRE1107, preconditionneur différent'])
legend('ilu','diag','none');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ A, matname] = orsirr_1;
n = size(A, 1);
b = (1:n)';
maxit = 2*n;
x0 = zeros(n,1);

[M11,M21] = ilu(A);

M12 = diag(diag(A));
M22 = eye(n);

M13 = eye(n);
M23 = M13';

subplot(2,1,2)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test avec m > iter == gmres

normRHS = norm(M21\(M11\b));
m = 30;
[x, flag, resrel, iter, resvec] = pGMRES(A, b, 0, tol, maxit, M11, M21, x0);
semilogy(resvec/normRHS);
hold on

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

normRHS = norm(M22\(M12\b));
m = 30;
[x, flag, resrel, iter, resvec] = pGMRES(A, b, 0, tol, maxit, M12, M22, x0);
semilogy(resvec/normRHS);
hold on
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

normRHS = norm(M23\(M13\b));
m = 30;
[x, flag, resrel, iter, resvec] = pGMRES(A, b, 0, tol, maxit, M13, M23, x0);
semilogy(resvec/normRHS);
hold on

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

semilogy(resvec/normRHS);
title(['ORSIRR1, preconditionneur différent'])
legend('ilu','diag','none');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
