function rhoE = RhoDE_2_RhoE(rhoDE, Ntau, Nc, displayFlag)
% Conversion d'une matrice densité énergie-délai (rhoDE)
% vers la matrice densité finale dans la base énergie pure (rhoE).
% Ntau : nombre de points de délai
% Nc   : nombre de canaux d'énergie
% displayFlag (optionnel) : true/false pour activer/désactiver l'affichage

    % Valeur par défaut : affichage activé si non précisé
    if nargin < 4
        displayFlag = true;
    end

    % Calcul du nombre d'échantillons effectif (NE) selon la convention "Sami"
    NE = Sampling_Sami(Nc, Ntau);

    % Décalage central pour indexation symétrique
    Ni = floor((Ntau-1)/2);

    % -- Construction d'une matrice densité augmentée (remplissage sur grand support)
    rhoErot_full = zeros(2*NE-1);
    % Insertion de rhoDE dans la fenêtre centrale de la grande matrice
    rhoErot_full(Ni+1:Ni+2*Nc, NE-(Ntau-1)/2 : NE+(Ntau-1)/2) = rhoDE;

    % -- Application d'un décalage cyclique (rotation) colonne par colonne
    rhoErot_fullfull = zeros(2*NE-1);
    for ic = 1:2*NE-1
        rhoErot_fullfull(:,ic) = circshift(rhoErot_full(:,ic), -ic + NE);
    end

    % --- Affichage intermédiaire de la densité "étalée" (si demandé)
    if displayFlag
        figure(819)
        imagesc(abs(rhoErot_fullfull))
        axis equal
        colormap(Carte_de_couleurs)
        title('|rhoErot\_fullfull| : densité énergie-délai étalée')
    end

    % -- Reconstruction finale de la matrice densité dans la base énergie
    rhoE = zeros(NE);
    for ir = 1:NE
        % V : vecteur d'indices pour chaque ligne (décalage le long de la diagonale)
        V = NE-ir+1 + (0:NE-1);
        rhoE(ir, :) = rhoErot_fullfull(2*ir-1, V);
    end

    % --- Affichages supplémentaires pour vérification (si demandé)
    if displayFlag
        figure(8181)
        subplot(1,3,1)
        imagesc(abs(rhoErot_full))
        axis equal
        title('|rhoErot\_full| : fenêtre insérée')

        subplot(1,3,2)
        imagesc(abs(rhoErot_fullfull))
        axis equal
        title('|rhoErot\_fullfull| : après rotation')

        subplot(1,3,3)
        imagesc(abs(rhoE))
        axis equal
        title('|rhoE| : matrice densité finale (base énergie)')
    end
end
