function [x, flag , relres, iter, resvec] = dqGMRES1(A, b,m, tol, maxit,M1,M2, x0)

%% Résolution de Ax =b par FOM ou GMRES
%% Paramètres
% A : la matrice du syst`eme que l’on cherche à résoudre
% b : le second membre de ce système
% x0 : le vecteur initial
% tol : le seuil demandé
% maxit : le nombre maximum d’itérations
%% Valeurs de retour
% x la solution
% flag qui indique si la méthode a convergé
% relres, la norme relative du résidu
% iter, le nombre d’it ́erations
% resvec, le vecteur des normes des résidus de chaque itération

%% Initialisation

flag = 1;
r =  M2\(M1\(b-A*x0));
beta = norm(r,2);
normR = beta;
normRHS = norm((M2\(M1\b)),2);
normB = norm(b,2);
relres = beta/normB;
v = [];
resvec = [];
resvec = [resvec beta];
v = [v r/beta];
j = 1;
x = x0;
convergence = false;
c = [];
s = [];
p = [];
gamma = beta;
gammanew = gamma;
%% Iteration
while (~convergence) && (j <= maxit)
    w = (M2\(M1\A))*v(:,j);
    %% Step1 : GS
    for i = max(1,j-m+1):j
        h(i,j) = v(:,i)'*w;
        w = w - h(i,j)*v(:,i);
    end
    h(j+1,j) = norm(w,2);
    v(:,j+1) = w/h(j+1,j);
    5
    
    %% Step 2 : facto qr
    for i = max(1,j-m):j-1
        h(i,j) = c(j)*h(i,j)+s(j)*h(i+1,j);
        h(i+1,j) = -s(j)*h(i,j)+c(j)*h(i+1,j);
    end
    
    %% Step 3 : Apply Qi to Hj
    c = [ c , h(j+1,j)/sqrt( h(j,j)^2 + h(j+1,j)^2)];
    s = [ s , h(j,j)/sqrt( h(j,j)^2 + h(j+1,j)^2)];
    gammanew = -s(j)*gamma;
    gamma = c(j)*gamma;
    h(j,j) = c(j)*h(j,j)+s(j)*h(j+1,j);
    %% Step 4
    aux = 0;
    for i = max(1,j-m):j-1
        aux = h(i,j)*p(i);
    end
    aux = 1/h(j,j)*(v(:,j)-aux);
    p = [ p , aux ];
    
    %% Test de convergence
    resvec = [resvec sqrt(j-m+1)*gammanew];
    relres = norm(resvec(end))/normRHS;
    convergence = (relres <= tol);
    j = j +1;
    gamma = gammanew;
end
iter = j;
if iter < maxit
    flag = 0;
end
end

