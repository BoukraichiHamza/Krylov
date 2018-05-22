function [x, flag , relres, iter, resvec] = dqGMRES1(A, b,m, tol, maxit,M1,M2, x0)

%% Résolution de Ax =b  dqGMRES incomplet
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
%% Iteration
while (~convergence) && (j <= maxit)
    w = (M2\(M1\A))*v(:,j);
    
    for i = max(1,j-m+1):j
        h(i,j) = v(:,i)'*w;
        w = w - h(i,j)*v(:,i);
    end
    h(j+1,j) = norm(w,2);
    v(:,j+1) = w/h(j+1,j);
    
    
    %% Facto qr
    [C,R] = qr(sparse(h(1:j+1,1:j)),(beta*eye(j+1,1)));
    y = R\C;
    
    x = x0 + v(:,1:j)*y;
    normR =  norm(C(end));
    resvec = [resvec normR];
    relres = normR/normRHS;
    convergence = (relres <= tol);
    j = j +1;
end
iter = j;
if iter < maxit
    flag = 0;
end
end

