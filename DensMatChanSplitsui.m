function [u] = DensMatChanSplitsui(A, W, u0, maxiter, tol, lmLambda)
% gaussNewtonRank1 Approxime A par u*u' sous ||W.*(A - u*u')||_F^2
%
% [u, history] = gaussNewtonRank1(A, W, u0, maxiter, tol, lmLambda)
%
% Entrées :
%   A        - matrice symétrique cible (n×n)
%   W        - matrice de poids (n×n), w_ij > 0
%   u0       - vecteur initial (n×1). Si vide, on prend le vecteur propre principal de A
%   maxiter  - nombre max. d’itérations (par défaut 100)
%   tol      - tolérance sur ||delta u||/||u|| pour la convergence (par défaut 1e-6)
%   lmLambda - paramètre de damping pour Levenberg–Marquardt (par défaut 1e-3)
%
% Sorties :
%   u        - vecteur tel que X = u*u' minimise ||W.*(A - X)||_F^2
%   history  - vecteur des valeurs de f(u) à chaque itération

    if nargin < 3 || isempty(u0),    u0 = [];           end
    if nargin < 4 || isempty(maxiter), maxiter = 100;    end
    if nargin < 5 || isempty(tol),     tol = 1e-6;        end
    if nargin < 6 || isempty(lmLambda), lmLambda = 1e-3;  end

    n = size(A,1);
    % Indices i<=j du triangle supérieur
    [I, J] = find(triu(true(n)));
    m = numel(I);

    % Initialisation de u
    if isempty(u0)
        B = (A + A')/2;
        [V, D] = eig(B);
        [~, idx] = max(diag(D));
        u = V(:, idx);
    else
        u = u0;
    end

    history = zeros(maxiter,1);

    for k = 1:maxiter
        % Calcul des résidus r (m×1)
        r = zeros(m,1);
        for t = 1:m
            i = I(t); j = J(t);
            r(t) = W(i,j) * (u(i)*u(j) - A(i,j));
        end

        % Objectif f = ||r||^2
        history(k) = r'*r;

        % Construction de la Jacobienne Jmat (m×n)
        Jmat = zeros(m,n);
        for t = 1:m
            i = I(t); j = J(t);
            Jmat(t,i) = W(i,j) * u(j);
            if i ~= j
                Jmat(t,j) = W(i,j) * u(i);
            end
        end

        % Équations normales avec damping LM
        JTJ  = Jmat' * Jmat + lmLambda * eye(n);
        grad = Jmat' * r;
        delta = - JTJ \ grad;

        % Mise à jour
        u_new = u + delta;

        % Critère d'arrêt
        if norm(delta) / norm(u) < tol
            u = u_new;
            history = history(1:k);
            break;
        end

        u = u_new;
    end
end
