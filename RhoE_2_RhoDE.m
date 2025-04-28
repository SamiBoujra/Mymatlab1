function rhoDE = RhoE_2_RhoDE(rhoE, Ntau, Nc, displayFlag)
    % Si le paramètre displayFlag n'est pas fourni, on l'active par défaut
    if nargin < 4
        displayFlag = true;
    end

    NE = size(rhoE, 1);
    rhoErot = zeros(2*NE-1);
    
    for ir = 1:NE
        V = NE-ir+1 + (0:NE-1);
        rhoErot(2*ir-1, V) = rhoE(ir, :);
    end

    % Réorganisation par décalage circulaire des colonnes
    for ic = 1:2*NE-1
        rhoErot(:, ic) = circshift(rhoErot(:, ic), ic-NE);
    end

    Ni = floor((Ntau-1)/2);
    rhoDE = rhoErot(Ni+1:Ni+2*Nc, NE-(Ntau-1)/2:NE+(Ntau-1)/2);

    % Affichage conditionnel
    if displayFlag
        figure(818)
        subplot(1,3,1)
        imagesc(abs(rhoE))
        axis equal
        title('abs(rhoE)')
        
        subplot(1,3,2)
        imagesc(abs(rhoErot))
        axis equal
        title('abs(rhoErot)')
        
        subplot(1,3,3)
        imagesc(abs(rhoDE))
        axis equal
        title('abs(rhoDE)')
    end
end