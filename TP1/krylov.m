function [x, flag , relres, iter, resvec] = krylov(A, b, x0, tol, maxit, type)

%% Résolution de Ax =b par FOM ou GMRES
%% Paramètres
% A : la matrice du syst`eme que l’on cherche à résoudre
% b : le second membre de ce système
% x0 : le vecteur initial
% tol : le seuil demandé
% maxit : le nombre maximum d’itérations
% type : 0 pour FOM, 1 pour GMRES
%% Valeurs de retour
% x la solution
% flag qui indique si la méthode a convergé
% relres, la norme relative du résidu
% iter, le nombre d’it ́erations
% resvec, le vecteur des normes des résidus de chaque itération

%% Initialisation

flag = 1;
r = b - A*x0;
beta = norm(r,2);
normb = norm(b,2);
relres = beta/normb;
v = [];
resvec = [];
resvec = [resvec beta];
v = [v r/beta];
j = 0;
%h = zeros(maxit+1,maxit);
x = x0;
%% Iteration
while (relres > tol) && (j <= maxit)
    j = j +1;
    w = A*v(:,j);
    
    for i = 1:j
        h(i,j) = v(:,i)'*w;
        w = w - h(i,j)*v(:,i);
    end
    h(j+1,j) = norm(w,2);
    v(:,j+1) = w/h(j+1,j);
    if (type == 1)
        
        y = h(1:j+1,1:j)\(beta*eye(j+1,1));
        x = x0 + v(:,1:j)*y;
    else
        y = h(1:j,1:j)\(beta*eye(j,1));
        x = x0 + v(:,1:j)*y;
    end
    
    r = b - A*x;
    normr = norm(r,2);
    resvec = [resvec normr];
    relres = normr/normb;
end
iter = j;
if iter < maxit
    flag = 0;
end
end

