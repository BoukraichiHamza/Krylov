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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test avec m > iter == gmres

subplot(2,2,1)
normRHS = norm(M21\(M11\b));
[x, flag, resrel, iter, resvec] = gmres(A, b, [], tol, maxit, M11, M21, x0);

fprintf(' flag : %5d \n', flag);
fprintf(' nbiter : %5d \n', iter);
fprintf(' résidu relatif: %3.1e \n', resrel);
semilogy(resvec/normRHS, 'x');
hold on;

m = 30;
[x, flag, resrel, iter, resvec] = pGMRES(A, b, 0, tol, maxit, M11, M21, x0);
fprintf(' flag : %5d \n', flag);
fprintf(' nbiter : %5d \n', iter);
fprintf(' résidu relatif: %3.1e \n', resrel);
semilogy(resvec/normRHS, 'o');
title(['mat1, tol = 1e-10, ilu precondtionner'],'FontSize',12);
legend('gmres', ['dqgmres(',num2str(m),')']);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,2)
normRHS = norm(M22\(M12\b));
[x, flag, resrel, iter, resvec] = gmres(A, b, [], tol, maxit, M12, M22, x0);

fprintf(' flag : %5d \n', flag);
fprintf(' nbiter : %5d \n', iter);
fprintf(' résidu relatif: %3.1e \n', resrel);
semilogy(resvec/normRHS, 'x');
hold on;

m = 30;
[x, flag, resrel, iter, resvec] = pGMRES(A, b, 0, tol, maxit, M12, M22, x0);
fprintf(' flag : %5d \n', flag);
fprintf(' nbiter : %5d \n', iter);
fprintf(' résidu relatif: %3.1e \n', resrel);
semilogy(resvec/normRHS, 'o');
title(['mat1, tol = 1e-10, diagonal'])
legend('gmres', ['pgmres(',num2str(m),')']);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,3)
normRHS = norm(M23\(M13\b));
[x, flag, resrel, iter, resvec] = gmres(A, b, [], tol, maxit, M13, M23, x0);

fprintf(' flag : %5d \n', flag);
fprintf(' nbiter : %5d \n', iter);
fprintf(' résidu relatif: %3.1e \n', resrel);
semilogy(resvec/normRHS, 'x');
hold on;

m = 30;
[x, flag, resrel, iter, resvec] = pGMRES(A, b, 0, tol, maxit, M13, M23, x0);
fprintf(' flag : %5d \n', flag);
fprintf(' nbiter : %5d \n', iter);
fprintf(' résidu relatif: %3.1e \n', resrel);
semilogy(resvec/normRHS, 'o');
title(['mat1, tol = 1e-10, none'])
legend('gmres', ['pgmres(',num2str(m),')']);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,4)
normRHS = norm(M24\(M14\b));
[x, flag, resrel, iter, resvec] = gmres(A, b, [], tol, maxit, M14, M24, x0);

fprintf(' flag : %5d \n', flag);
fprintf(' nbiter : %5d \n', iter);
fprintf(' résidu relatif: %3.1e \n', resrel);
semilogy(resvec/normRHS, 'x');
hold on;

m = 30;
[x, flag, resrel, iter, resvec] = pGMRES(A, b, 0, tol, maxit, M14, M24, x0);
fprintf(' flag : %5d \n', flag);
fprintf(' nbiter : %5d \n', iter);
fprintf(' résidu relatif: %3.1e \n', resrel);
semilogy(resvec/normRHS, 'o');
title(['mat1, tol = 1e-10, chol'])
legend('gmres', ['pgmres(',num2str(m),')']);