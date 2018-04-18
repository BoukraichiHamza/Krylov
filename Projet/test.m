close all
clear all
mat1 = load('mat1.mat');
pde225_5e1 = load('pde225_5e-1.mat');
hydcar20 = load('hydcar20.mat');

Mat1 = mat1.A;
bmat1 = [1:size(Mat1,1)]';

Mat2 = pde225_5e1.A;
bmat2 = [1:size(Mat2,1)]';

x0 = zeros(size(Mat2,1),1);

M1 = diag(diag(sqrt(Mat2)));
M2 = M1;

% [x1, flag1, relres1, iter1, resvec1] = pGMRES(Mat1, bmat1, 0, 1e-6,500,M1,M2,x0);
% [x11, flag11, relres11, iter11, resvec11] = gmres(Mat1, bmat1, [], 1e-6,500,M1,M2,x0);
% 
% 
% figure(1)
% semilogy(resvec1/norm(bmat1,2),'o')
% hold on
% semilogy(resvec11/norm(bmat1,2),'*')
%i = norm(resvec11 - resvec1',Inf)
m = 20;
[x1, flag1, relres1, iter1, resvec1] = dqGMRES2(Mat2, bmat2, m, 1e-6,500,M1,M2,x0);
[x11, flag11, relres11, iter11, resvec11] = gmres(Mat2, bmat2, m, 1e-6,500,M1,M2,x0);


figure(1)
semilogy(resvec1/norm(bmat1,2),'o')
hold on
semilogy(resvec11/norm(bmat1,2),'*')