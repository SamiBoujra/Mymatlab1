function rhoDE = RhoE_2_RhoDE(rhoE, Ntau, Nc, displayFlag)
    % Conversion de la matrice densité dans la base énergie (rhoE)
    % vers la base énergie-délai (rhoDE).
    % Ntau : nombre de points de délai
    % Nc   : nombre de canaux d'énergie
    % displayFlag : (optionnel) true pour afficher les étapes

    % Si le paramètre displayFlag n'est pas fourni, on l'active par défaut
    if nargin < 4
        displayFlag = true;
    end

    % Taille de la base énergie (nombre de points)
    NE = size(rhoE, 1);

    % Initialisation d'une grande matrice temporaire pour la "rotation"
    rhoErot = zeros(2*NE-1);

    % --- Remplissage de la grande matrice par lignes décalées ---
    % Pour chaque ligne de rhoE, on la place sur la ligne (2*ir-1)
    % et on la décale horizontalement selon l'indice ir
    for ir = 1:NE
        % V : indices colonnes de la fenêtre associée à la ligne ir
        V = NE-ir+1 + (0:NE-1);
        rhoErot(2*ir-1, V) = rhoE(ir, :);
    end

    % --- Décalage circulaire de chaque colonne ---
    % Permet d'aligner correctement les contributions énergie-délai
    for ic = 1:2*NE-1
        % circshift vertical des colonnes (rotation)
        rhoErot(:, ic) = circshift(rhoErot(:, ic), ic-NE);
    end

    % --- Extraction de la fenêtre centrale correspondant à (Nc,Ntau) ---
    % On sélectionne la sous-matrice qui contient les valeurs utiles pour rhoDE
    Ni = floor((Ntau-1)/2);
    rhoDE = rhoErot(Ni+1:Ni+2*Nc, NE-(Ntau-1)/2 : NE+(Ntau-1)/2);

    % --- Affichage conditionnel pour vérification ---
    if displayFlag
        figure(818)
        subplot(1,3,1)
        imagesc(abs(rhoE))
        axis equal
        title('abs(rhoE)') % Matrice densité énergie initiale

        subplot(1,3,2)
        imagesc(abs(rhoErot))
        axis equal
        title('abs(rhoErot)') % Matrice densité "étalée" après rotation

        subplot(1,3,3)
        imagesc(abs(rhoDE))
        axis equal
        title('abs(rhoDE)') % Matrice densité énergie-délai finale
    end
end
