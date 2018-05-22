function [x, flag , relres, iter, resvec] = dqGMRES1(A, b,m, tol, maxit,M1,M2, x0)

%% Résolution de Ax =b par dqGMRES
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
j = 0;
x = x0;
convergence = false;
c = zeros(m,1);
s = zeros(m,1);
p = [];
gamma = beta;
gammanew = gamma;
h = zeros(maxit+1,maxit);
%% Iteration
while (~convergence) && (j < maxit)
    j = j + 1;
    w = (M2\(M1\A))*v(:,j);
    %% Step1 : GS
    for i = max(1,j-m+1):j
        h(i,j) = v(:,i)'*w;
        w = w - h(i,j)*v(:,i);
    end
    h(j+1,j) = norm(w,2);
    v(:,j+1) = w/h(j+1,j);
    
    %% Step 2 : facto qr
    for i = max(1,j-m):j-1
        newhij = c(i)*h(i,j)+s(i)*h(i+1,j);
        newhi1j = -s(i)*h(i,j)+c(i)*h(i+1,j);
        h(i,j) = newhij;
        h(i+1,j) = newhi1j;
    end

    %% Step 3 : Apply Qi to Hj
    c(j) =  h(j,j)/sqrt( h(j,j)^2 + h(j+1,j)^2);
    s(j) =  h(j+1,j)/sqrt( h(j,j)^2 + h(j+1,j)^2);
    gammanew = -s(j)*gamma;
    gamma = c(j)*gamma;
    h(j,j) = c(j)*h(j,j)+s(j)*h(j+1,j);
    h(j+1,j) = 0;
    %% Step 4
    p(:,j) = v(:,j)/h(j,j);
    for i = max(1,j-m):j-1
        p(:,j) = p(:,j)-h(i,j)*p(:,i)/h(j,j);
    end
    x =x+gamma*p(:,j);
    
    %% Test de convergence
    resvec = [resvec abs(gammanew)];
    relres = resvec(end)/normRHS;
    convergence = (relres <= tol);
    gamma = gammanew;
end
iter = j;
if iter < maxit
    flag = 0;
end
end

