close all;
clc;
clear;

%% ===========================
% 1. Initialisation du pool parallèle
% ===========================
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool; % Lancement d'un pool avec le nombre de workers par défaut
end

%% ===========================
% 2. Chargement des constantes physiques
% ===========================
% Les variables (eV, h, h_bar, fs, c, epsilon0, TW, etc.) sont chargées depuis constants.mat
load constants.mat

%% ===========================
% 3. Axes d'énergies et de délais
% ===========================
Ntau = 150;             % Nombre de points sur l'axe du délai
dtau = 0.01 * fs;       % Pas de l'axe du délai (fs)
N = 256;                % Nombre de points sur l'axe énergie/cinétique
E0 = 10;                % Énergie minimale (eV)
DeltaE = 80;            % Largeur de la fenêtre énergétique (eV)
dt = h/(DeltaE*eV);     % Pas temporel associé via la relation d'incertitude
dE = DeltaE/N;          % Pas énergétique
KE = E0 + (0:N-1)*dE;   % Axe des énergies cinétiques

% Trouve Zmax pour que Nc ≤ N (contrainte de troncature)
DeltaEf = 80;
zCandidate = 0;
while true
    tNc = round(2*(zCandidate+1)*DeltaEf*(eV/h)*Ntau*dtau);
    if tNc <= N
        zCandidate = zCandidate + 1;
    else
        break;
    end
end
Z = max(0, zCandidate-1);

Nc = round(2 * Z * DeltaEf * (eV/h) * Ntau * dtau); % Nombre de canaux effectif
Nf = Nc + ceil((Ntau-1)/2);      % Décalage associé au delay
NE = 2*Nc;                       % Taille totale pour les matrices (double pour symétrisation)

%% ===========================
% 4. Paramètres des harmoniques (XUV)
% ===========================
HH = 39:2:47;        % Ordres harmoniques
NHH = numel(HH);     % Nombre d'harmoniques

%% ===========================
% 5. Décorrélation/Décohérence diagonale (simule le bruit laser etc.)
% ===========================
NDtP = 250;                  % Nb points sur l’axe du délai
nDtP = 0:NDtP-1;             % Axe du délai
DiagDecoh_qqp = zeros(2*NDtP-1, NHH, NHH);
for iq = 1:NHH
    for iqp = 1:NHH
        FWHM = round(NDtP*2/(HH(iq)+HH(iqp))); % Largeur à mi-hauteur
        G = PulseShape('gauss', [2*NDtP-1,FWHM], 1); % Fonction gaussienne centrée
        DiagDecoh_qqp(:,iq,iqp) = G./sum(G);  % Normalisation
    end
end

% Axes énergie centrale pour la densité d'état
dEc = DeltaEf/Nc;
E0f = 10;     % Centre énergétique (eV)
E0c = E0f;
nc = 0:Nc-1;
Ec = E0c + nc*dEc;   % Axe énergie centrale

%% ===========================
% 6. Calcul du champ IR (modulateur de phase)
% ===========================
E0_IR = au2SI(0.057, 'energy') / eV;    % Energie photon IR (en eV)
T_IR = h / (E0_IR*eV);                  % Période de l’IR
omega = 2*pi/T_IR;                      % Fréquence angulaire

shapeI = PulseShape('gauss', [N, round(3*fs*sqrt(2)/dt)], 1);  % Profil du pulse IR
t = (0:N-1)*dt;    % Axe temporel
A_IR0 = -sin(omega*t + pi/2).' .* shapeI(:) / omega;   % Potentiel vecteur IR

%% ===========================
% 7. Final sampling et réglages
% ===========================
Nzp = 0;                 % Zero-padding éventuel
Nt_P = NDtP + Nzp;       % Taille totale avec padding
Nmodes = 6;              % Modes spectraux à conserver
HH = 39:2:45;            % Ordres harmoniques utilisés pour la simulation
I_HH = [0.25 0.5 0.75 1 0.1];  % Intensité relative de chaque harmonique
NHH = length(HH);
Phi = [-10 +0.25 -0.4 -4.5 ];  % Phases relatives
dE = DeltaE / Nt_P;
E = E0 + (0:Nt_P-1)*dE;
t = (-Nt_P/2:Nt_P/2-1)*dt;
shape = PulseShape('sin2', [NDtP NDtP/2], 1);
Ip = GiveIp('Ne');       % Potentiel d’ionisation du Néon

%% ===========================
% 8. Calcul de la densité de matrice d’état (rho)
% ===========================
reso = 1;
if reso
    % Réinitialisation de la décohérence diagonale pour ce cas
    DiagDecoh_qqp = zeros(2*NDtP-1, NHH, NHH);
    for iq = 1:NHH
        for iqp = 1:NHH
            FWHMqqp = round(NDtP*2/(iq+iqp));
            DiagDecoh_qqp(:,iq,iqp) = PulseShape('gauss', [2*NDtP-1 FWHMqqp], 1);
            DiagDecoh_qqp(:,iq,iqp) = DiagDecoh_qqp(:,iq,iqp)./sum(DiagDecoh_qqp(:,iq,iqp));
        end
    end
    rhoDtP = zeros(NDtP);
    for iq = 1:NHH
        for iqp = 1:NHH
            rhot_qqp = sqrt(I_HH(iq))*sqrt(I_HH(iqp))*exp( ...
                -1i*2*pi*(E0_IR*HH(iq)-E0-Ip)*eV/h*nDtP(:)*dt - ...
                1i*(Phi(iq)-Phi(iqp)) + ...
                1i*2*pi*(E0_IR*HH(iqp)-E0-Ip)*eV/h*nDtP*dt ...
                );
            rhot_qqp = AddDiagDecoh(rhot_qqp,DiagDecoh_qqp(:,iq,iqp));
            rhoDtP = rhoDtP + rhot_qqp;
        end
    end
    rhoDtP = rhoDtP.*(shape(:)*shape);
    rhoDtP = rhoDtP/trace(rhoDtP);
    rhot = zeros(Nt_P);
    rhot(Nzp/2+1:Nzp/2+NDtP, Nzp/2+1:Nzp/2+NDtP) = rhoDtP;
    rhoE = SwitchStateRep(rhot, 'Rhot2RhoE', 1);  % Passage vers la base énergie
    Purity(rhot);
    sigma = eig(rhot);
    figure(1414); semilogy(abs(sigma), 'ro-');
    [V,D] = eigs(rhot, Nmodes);
    P = V*sqrt(D); RhoEig = P*P';
    nu = Purity(RhoEig);
    rhoEEig = SwitchStateRep(RhoEig, 'Rhot2RhoE', 1);
    figure(1515); imagesc((abs(rhoEEig))); colormap(Carte_de_couleurs)
else
    P = zeros(NDtP,1);
    for iHH = 1:NHH
        P = P + sqrt(I_HH(iHH))*exp(-1i*2*pi*(E0_IR*HH(iHH)-E0-Ip)*eV/h*nDtP(:)*dt);
    end
    P = P.*shape(:);
    P = vertcat(zeros(Nzp/2,1),P,zeros(Nzp/2,1));
end

rhoE = SwitchStateRep(rhot, 'Rhot2RhoE', 1);
rhoDE = RhoE_2_RhoDE(rhoE, Ntau, Nc);
rhoE2 = RhoDE_2_RhoE(rhoDE, Ntau, Nc);

%% ===========================
% 9. Préparation des opérateurs de modulation et des matrices physiques
% ===========================
E0_IR = 1.55;  % Energie centrale du photon laser IR (eV)
T_IR = h/(E0_IR*eV);
omega_IR = 2*pi/T_IR;
t = (-N/2:N/2-1)*dt;
p = KE2p(KE);
Intensity = 8;   % Intensité IR en TW/cm2
IntensityW = Intensity*1e4*TW; % En W/m2
FWHM_IR_intensity = 3*fs;
FWHM_IRfield = FWHM_IR_intensity*sqrt(2);

shape = PulseShape('gauss', [N, round(FWHM_IRfield/dt)], 1);
A_IR0 = -sin(omega_IR*t(:)+pi/2).*shape(:)/omega_IR;
A_IR = A_IR0*sqrt(2*IntensityW/(c*epsilon0)); % Champ IR effectif

[Phi_mod, intA, intA2] = Phase_modulator(A_IR, dt, p);
Gt = exp(1i*Phi_mod);

% Transformation temporelle → énergétique
Gtild = (ifft(Gt, [], 2));
Gtild = fftshift(Gtild.', 1);
Ef = Ec;

%% ===========================
% 10. Construction des matrices de convolution / opérateurs de reconstruction
% ===========================
Ctot = zeros(N, 2*NE, Ntau);
parfor ipf = 1:N
    g_vec = Gtild(:, ipf);          % Vecteur col énergie
    C = g_vec * g_vec.';            % Produit extérieur
    Crot = fftshift(RhoE_2_RhoDE(C, Ntau, Nc, false), 2); % Passage DE, rotation
    Cext = [Crot; zeros(NE, Ntau)]; % Extension (padding)
    Csh = fftshift(circshift(Cext, [ipf-1, 0]));
    Ctot(ipf, :, :) = Csh;
end

%% ===========================
% 11. Construction du sinogramme
% ===========================
Stild = zeros(N, Ntau);
for iDE = 1:Ntau
    vRhoDE = rhoDE(:,iDE);
    % Partitionne Ctot en 2 pour chaque retard
    M1 = squeeze(Ctot(:,1:NE/2,iDE));
    M2 = squeeze(Ctot(:,NE/2+NE+1:end,iDE));
    % Concatène, FFT-shift, puis applique la densité
    M = fftshift([M1,M2],2);
    M_new = M;
    Stild(:,iDE) = M_new*vRhoDE;
end
Stild(isnan(Stild)) = 0; % Sécurité NaN

%% ===========================
% 12. Ajout de bruit (optionnel)
% ===========================
S = fft(Stild, [], 2);
noise_level = 0e-6;
noise = (randn(size(S)) + 1i*randn(size(S)))/sqrt(2);
modS = abs(S);
phiS = angle(S);
modS_noisy = modS + noise_level * randn(size(S));
S_noisy = max(modS_noisy, 0) .* exp(1i * phiS);
Stild2 = ifft(S_noisy, [], 2);

figure;
subplot(2,1,1); imagesc(abs(S)); colorbar; title('Ideal S');
subplot(2,1,2); imagesc(abs(Stild2)); colorbar; title('Noisy S');

figure;
subplot(2,1,1); imagesc(angle(S_noisy)); colorbar; title('Noisy S (angle)');
subplot(2,1,2); imagesc(angle(S)); colorbar; title('Ideal S (angle)');

%% ===========================
% 13. Préparation des opérateurs pour la reconstruction (gradients, Laplacien)
% ===========================
NE = 348; e = ones(NE,1);
D1 = spdiags([-e e], [0 1], NE, NE);
D2 = spdiags([e -2*e e], [-1 0 1], NE, NE);
W2 = D1'*D1; L2 = D2'*D2;

g1_list = logspace(-6,5,20);
g2_list = logspace(-6,5,20);
g3_list = logspace(-4,4,8);
iter_max = 100;

%% ===========================
% 14. Reconstruction module + phase (Tikhonov + TV)
% ===========================
E = ones(NE,1);
D1 = spdiags([-E E], [0 1], NE-1, NE);
D2 = spdiags([E -2*E E], [-1 0 1], NE, NE);

alpha1 = 1e-4; alpha2 = 1e-4;
Reg_block = alpha1 * blkdiag(D1'*D1, D1'*D1) + alpha2 * blkdiag(D2'*D2, D2'*D2);

rhoDE_modulus = zeros(NE, Ntau);
phi_init = zeros(NE, Ntau);

parfor iDE = 1:Ntau
    % a) Récupère la matrice d’opérateur M pour le delay courant
    M1 = squeeze(Ctot(:,1:NE/2, iDE));
    M2 = squeeze(Ctot(:,NE/2+NE+1:end, iDE));
    M = fftshift([M1, M2], 2);
    d = Stild2(:, iDE);

    % b) Nettoyage SVD (troncature basse)
    [U,S2,V] = svd(M,'econ'); s = diag(S2);
    r = max(find(s>=0.01*s(1),1,'last'),10);
    Mclean = U(:,1:r)*S2(1:r,1:r)*V(:,1:r)';

    % c) Résolution linéaire complexe (réel/imaginaire empilé)
    A_stack = [real(Mclean), -imag(Mclean); imag(Mclean), real(Mclean)];
    b_stack = [real(d); imag(d)];
    x_sol = (A_stack'*A_stack + Reg_block) \ (A_stack'*b_stack);
    x_cplx = x_sol(1:NE) + 1i*x_sol(NE+1:end);

    % d) Filtrage-médiane du module (suppression du background)
    abs_m = abs(x_cplx);
    win = max(5, round(NE/10));
    baseline = medfilt1(abs_m, win);
    m_ref = abs_m - baseline; m_ref(m_ref < 0) = 0;
    thr = 0.2 * max(m_ref); m_ref(m_ref < thr) = 0;
    rhoDE_modulus(:,iDE) = abs(m_ref);
    phi_init(:,iDE) = angle(x_cplx);
end

rhoDE_phase = phi_init;

%% ===========================
% 15. Passage densité état finale
% ===========================
rhoDE_weighted = rhoDE_modulus .* exp(1i * phi_init);
rhoE_weighted = RhoDE_2_RhoE(rhoDE_weighted, Ntau, Nc);

% Visualisation
figure(138);
subplot(1,2,1); imagesc(abs(rhoE_weighted)); axis image; colorbar; title('|ρ| final');
subplot(1,2,2); imagesc(angle(rhoE_weighted)); axis image; colorbar; title('∠ρ low-rank');

figure(139)
subplot(1,2,1); imagesc(abs(rhoE)); axis image; colorbar; title('|ρ| initial');
subplot(1,2,2); imagesc(angle(rhoE)); axis image; colorbar; title('∠ρ initial');

%% ===========================
% 16. Bootstrap par rééchantillonnage des résidus
% ===========================
Nboot = 20;
rhoDE_modulus_boot = zeros(NE, Ntau, Nboot);
phi_init_boot = zeros(NE, Ntau, Nboot);

E = ones(NE,1);
D1 = spdiags([-E E], [0 1], NE-1, NE);
D2 = spdiags([E -2*E E], [-1 0 1], NE, NE);
alpha1 = 1e-4; alpha2 = 1e-4;
Reg_block = alpha1 * blkdiag(D1'*D1, D1'*D1) + alpha2 * blkdiag(D2'*D2, D2'*D2);

phi_init = zeros(NE, Ntau);
rhoDE_modulus = zeros(NE, Ntau);
d_hat_store = zeros(size(Stild2));

% Reconstruction initiale, stockage du signal reconstruit d_hat
parfor iDE = 1:Ntau
    M1 = squeeze(Ctot(:,1:NE/2, iDE));
    M2 = squeeze(Ctot(:,NE/2+NE+1:end, iDE));
    M = fftshift([M1, M2], 2);
    d = Stild2(:, iDE);

    % SVD tronquée pour la stabilité numérique
    [U,S2,V] = svd(M,'econ'); s = diag(S2);
    r = max(find(s>=0.01*s(1),1,'last'),10);
    Mclean = U(:,1:r)*S2(1:r,1:r)*V(:,1:r)';

    % Résolution linéaire complexe (réel/imaginaire empilé)
    A_stack = [real(Mclean), -imag(Mclean); imag(Mclean), real(Mclean)];
    b_stack = [real(d); imag(d)];
    x_sol = (A_stack'*A_stack + Reg_block) \ (A_stack'*b_stack);
    x_cplx = x_sol(1:NE) + 1i*x_sol(NE+1:end);

    % Filtrage médian (suppression du bruit de fond)
    abs_m = abs(x_cplx);
    win = max(5, round(NE/10));
    baseline = medfilt1(abs_m, win);
    m_ref = abs_m - baseline; m_ref(m_ref < 0) = 0;
    thr = 0.2 * max(m_ref); m_ref(m_ref < thr) = 0;

    % Stockage module/phase bruts et signal reconstruit
    rhoDE_modulus(:,iDE) = abs(x_cplx);
    phi_init(:,iDE) = angle(x_cplx);
    rhoDE_raw(:,iDE) = x_cplx;
    d_hat_store(:,iDE) = M * x_cplx;
end

%% === 2) Bootstrap sur les résidus ===
rhoDE_boot = zeros(NE, Ntau, Nboot);    % sinogrammes bootstrap
rhoE_boot = zeros(250, 250, Nboot);       % matrices densité EE bootstrap

for b = 1:Nboot
    parfor iDE = 1:Ntau
        % Même pipeline que ci-dessus mais sur des données "perturbées"
        M1 = squeeze(Ctot(:,1:NE/2, iDE));
        M2 = squeeze(Ctot(:,NE/2+NE+1:end,iDE));
        M = fftshift([M1, M2],2);

        d0 = Stild2(:, iDE);
        d_hat = d_hat_store(:, iDE);
        % Bootstrap : rééchantillonnage des résidus
        res = d0 - d_hat;
        d_star = d_hat + res(randsample(numel(res), numel(res), true));

        % SVD tronquée
        [U,S2,V] = svd(M,'econ'); s=diag(S2);
        r = max(find(s>=0.01*s(1),1,'last'),10);
        Mclean = U(:,1:r)*S2(1:r,1:r)*V(:,1:r)';

        % Résolution complexe
        A_stack = [ real(Mclean), -imag(Mclean); imag(Mclean), real(Mclean) ];
        b_stack = [ real(d_star); imag(d_star) ];
        x_sol = (A_stack'*A_stack + Reg_block) \ (A_stack'*b_stack);
        x_cplx = x_sol(1:NE) + 1i*x_sol(NE+1:end);

        % Nettoyage module (filtrage médian et seuillage)
        abs_m = abs(x_cplx);
        win = max(5,round(NE/10));
        baseline = medfilt1(abs_m,win);
        m_ref = abs_m - baseline; m_ref(m_ref<0)=0;
        thr = 0.2*max(m_ref); m_ref(m_ref<thr)=0;

        % Stocke sinogramme (on garde la phase brute)
        phi_ref = angle(x_cplx);
        rhoDE_boot(:,iDE,b) = m_ref .* exp(1i*phi_ref);
    end

    % === 2) Passage DE→E, projection semi-définie positive
    sinogram = rhoDE_boot(:,:,b);
    rhoE_rec = RhoDE_2_RhoE(sinogram, Ntau, Nc);
    rhoE_rec = (rhoE_rec + rhoE_rec')/2; % symétrisation hermitienne

    % Décalage spectral minimal pour garantir PSD
    [Vb,Db] = eig(rhoE_rec);
    db = diag(Db);
    tol = 1e-12;
    delta = max(0, -min(db) + tol);
    rhoE_pd = Vb * (Db + delta*eye(size(Db))) * Vb';

    % Renormalisation trace=1
    rhoE_pd = rhoE_pd / trace(rhoE_pd);
    rhoE_boot(:,:,b) = rhoE_pd;

    % Contrôle valeurs propres
    [~,D] = eig(rhoE_boot(:,:,b));
    eigvals = diag(D);
    tol_eig = -1e-12;
    if all(eigvals > tol_eig)
        disp(['✅ Bootstrap ', num2str(b), '/', num2str(Nboot), ' PSD OK.']);
    else
        warning('⚠️ Certaines valeurs propres sont négatives : min(λ) = %g', min(eigvals));
    end
end

%% === 3) Post-traitement Bootstrap : Moyenne, Écart-type, Visualisation ===
rho_mod_mean = mean(rhoE_boot, 3);
rho_mod_std  = std(rhoE_boot,  0, 3);

figure(130);
subplot(1,2,1); imagesc(abs(rho_mod_mean)); axis image; colorbar; title('|ρ| : moyenne bootstrap');
subplot(1,2,2); imagesc(angle(rho_mod_mean)); axis image; colorbar; title('Arg(ρ) : moyenne bootstrap');

figure(131);
subplot(1,2,1); imagesc(abs(rho_mod_std)); axis image; colorbar; title('|ρ| : écart-type bootstrap');
subplot(1,2,2); imagesc(angle(rho_mod_std)); axis image; colorbar; title('Arg(ρ) : écart-type bootstrap');

%% === FIN DU SCRIPT PRINCIPAL ===
