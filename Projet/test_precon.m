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
figure
subplot(2,2,1)

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
title(['mat1, preconditionneur différent'])
legend('ilu','diag','none','ichol');