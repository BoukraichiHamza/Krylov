close all
clear all
load('mat1.mat');
% pde225_5e1 = load('pde225_5e-1.mat');
% hydcar20 = load('hydcar20.mat');

% [ Mat1 , matname] = gre_1107;
% bmat1 = [1:size(Mat1,1)]';
normRHS = norm(b);
n = size(A,1);
x0 = zeros(n,1);
% M1 = diag(diag(sqrt(Mat1)));
% M2 = M1;

[M1,M2] = ilu(A);

 [x1, flag1, relres1, iter1, resvec1] = restartedGMRES(A, b, 30, 1e-10,2*n,M1,M2,x0);
 [x11, flag11, relres11, iter11, resvec11] = pGMRES(A, b,0, 1e-10,2*n,M1,M2,x0);
 [x12, flag12, relres12, iter12, resvec12] = dqGMRES2(A, b,30, 1e-10,2*n,M1,M2,x0);

 figure(1)
 semilogy(resvec1/normRHS)
 hold on
 semilogy(resvec11/normRHS)
 hold on
 semilogy(resvec12/normRHS)
 hold on
 title('Comparaison des 3 méthodes')
 legend('restartedGMRES','pGMRES','dqGMRES')
% [x1, flag1, relres1, iter1, resvec1] = dqGMRES2(Mat1, bmat1, m, 1e-6,max_it,M1,M2,x0);
% [x11, flag11, relres11, iter11, resvec11] = gmres(Mat1, bmat1, [], 1e-6,max_it,M1,M2,x0);
% 
% 
% figure(2)
% semilogy(resvec1/norm(bmat1,2),'o')
% hold on
% semilogy(resvec11/norm(bmat1,2),'*')