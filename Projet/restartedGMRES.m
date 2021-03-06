function [x, flag , relres, iter, resvec] = restartedGMRES(A, b,m, tol, maxit,M1,M2, x0)

%% Résolution de Ax =b par restartedGMRES
%% Paramètres
% A : la matrice du syst`eme que l’on cherche à résoudre
% b : le second membre de ce système
% x0 : le vecteur initial
% tol : le seuil demandé
% maxit : le nombre maximum d’itérations
% M1,M2 : composante du conditionneur
% x0 : itéré de départ
%% Valeurs de retour
% x : la solution
% flag : qui indique si la méthode a convergé
% relres : la norme relative du résidu
% iter: nombre d'itération
% resvec : vecteur des normes des résidus

convergence = false;
iter = 1;
x =x0;
resvec = [norm(M2\(M1\(b-A*x0)))];
while(~convergence) && (iter <= maxit)
    [x, flag , relres, iteraux, resvecaux] = pGMRES(A, b,0, tol, m,M1,M2, x);
    convergence = ( flag == 0);
    iter = iter + iteraux;
    resvec = [ resvec resvecaux(2:end)];
end

