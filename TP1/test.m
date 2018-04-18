load('mat1.mat');
n = size(A,1);
b = [1:n]';
x0 = zeros(n,1);
maxit = n;
tol = 1e-6;
[x, flag, relres, iter, resvec] = krylov(A, b, x0, tol, maxit, 1);
figure(1)
semilogy(resvec/norm(b,2))
hold on
[xm,flagm,relersm,iterm,resvecm] = gmres(A,b,n,tol,maxit);
semilogy(resvecm/norm(b,2))
flag
flagm
norm(x - xm,2)
iter
iterm

load('hydcar20.mat');
n = size(A,1);
b = [1:n]';
x0 = zeros(n,1);
maxit = n;
tol = 1e-6;
[x, flag, relres, iter, resvec] = krylov(A, b, x0, tol, maxit, 0);

figure(2)
semilogy(resvec/norm(b,2))



