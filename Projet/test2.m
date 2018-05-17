clear all;
close all;
load mat1.mat;

tol = 1e-10;
n = size(A, 1);
maxit = 2*n;
x0 = zeros(n,1);

[M1,M2] = ilu(A);
%M1 = diag(diag(A));
%M1 = eye(n);
%M2 = M1';

normRHS = norm(M2\(M1\b));
figure

% test avec m > iter == gmres
[x, flag, resrel, iter, resvec] = gmres(A, b, [], tol, maxit, M1, M2, x0);

fprintf(' flag : %5d \n', flag);
fprintf(' nbiter : %5d \n', iter);
fprintf(' résidu relatif: %3.1e \n', resrel);
semilogy(resvec/normRHS, 'x');
hold on;

m = 30;
[x, flag, resrel, iter, resvec] = dqGMRES2(A, b, m, tol, maxit, M1, M2, x0);
fprintf(' flag : %5d \n', flag);
fprintf(' nbiter : %5d \n', iter);
fprintf(' résidu relatif: %3.1e \n', resrel);
semilogy(resvec/normRHS, 'o');
title(['mat1, tol = 1e-10, ilu precondtionner'],'FontSize',12);
legend('gmres', ['dqgmres(',num2str(m),')']);
