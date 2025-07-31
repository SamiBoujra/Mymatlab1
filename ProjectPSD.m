function rho_psd = ProjectPSD(rho, tol)
% Projet une matrice sur l'ensemble des matrices semi-définies positives
% - rho : matrice de densité potentiellement non-PSD
% - tol : seuil de tolérance (défaut 1e-12)
%
% Sortie : rho_psd (PSD, Hermitienne, trace=1)

if nargin < 2
    tol = 1e-12;
end

% Forcer hermitienne
rho = (rho + rho') / 2;

% Diagonalisation
[V, D] = eig(rho);

% Prendre seulement la partie réelle des valeurs propres (partie imaginaire = bruit num)
d = real(diag(D));

% Remplace toute valeur négative ou trop petite par zéro
d(d < tol) = 0;

% Reconstruire matrice
rho_psd = V * diag(d) * V';

% Forcer hermitienne après reconstruction (précaution)
rho_psd = (rho_psd + rho_psd')/2;

% Renormalise la trace
rho_psd = rho_psd / trace(rho_psd);

end