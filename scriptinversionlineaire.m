close all;
clc;
clear;
%% Initialisation du pool parallèle
poolobj = gcp('nocreate');
if isempty(poolobj)
parpool; % lance un pool avec le nombre de workers par défaut
end
%% Chargement des constantes
load constants.mat
load constants.mat % eV, h, h_bar, fs, c, epsilon0, TW, etc.
% Delay axis
Ntau = 150;
dtau = 0.01 * fs;
% KE axis
N = 256;
E0 = 10; DeltaE = 80;
dt = h/(DeltaE*eV);
dE = DeltaE/N;
KE = E0 + (0:N-1)*dE;
% Trouve Zmax pour que Nc ≤ N
DeltaEf = 80;
zCandidate = 0;
while true
tNc = round(2*(zCandidate+1)*DeltaEf*(eV/h)*Ntau*dtau);
if tNc <= N, zCandidate = zCandidate+1;
else break; end
end
Z = max(0, zCandidate-1);
% Delay axis
Ntau = 150;
dtau = 0.01 * fs;
Nc = round(2 * Z * DeltaEf * (eV/h) * Ntau * dtau);
Nf =Nc+ceil((Ntau-1)/2);
NE =2*Nc;
% Harmoniques
HH = 39:2:47;
NHH = numel(HH);
% Décohérence
NDtP = 250;
nDtP = 0:NDtP-1;
DiagDecoh_qqp = zeros(2*NDtP-1, NHH, NHH);
for iq = 1:NHH
for iqp = 1:NHH
FWHM = round(NDtP*2/(HH(iq)+HH(iqp)));
G = PulseShape('gauss',[2*NDtP-1,FWHM],1);
DiagDecoh_qqp(:,iq,iqp) = G./sum(G);
end
end
dEc = DeltaEf/Nc;
E0f = 10;%eV
E0c = E0f;
nc = 0:Nc-1;
Ec = E0c + nc*dEc;
% Pré-calcule A_IR0 pour Phase_modulator
E0_IR = au2SI(0.057,'energy')/eV;
T_IR = h/(E0_IR*eV);
omega = 2*pi/T_IR;
shapeI = PulseShape('gauss',[N, round(3*fs*sqrt(2)/dt)],1);
t = (0:N-1)*dt;
A_IR0 = -sin(omega*t + pi/2).' .* shapeI(:) / omega;
%% Final sampling
Nzp = 0;%nombre pair
%Nzp = 156;%nombre pair
Nt_P = NDtP + Nzp;
Nmodes = 6;
HH = 39:2:45;
I_HH = [0.25 0.5 0.75 1 0.1];
NHH = length(HH);
Phi = [-10 +0.25 -0.4 -4.5 ];
dE = DeltaE/Nt_P;
E = E0 + (0:Nt_P-1)*dE;
t = (-Nt_P/2:Nt_P/2-1)*dt;
shape = PulseShape('sin2',[NDtP NDtP/2],1);
Ip = GiveIp('Ne');
%% Resolution
reso = 1;
if reso
DiagDecoh_qqp = zeros(2*NDtP-1,NHH,NHH);
for iq = 1:NHH
for iqp = 1:NHH
FWHMqqp = round(NDtP*2/(iq+iqp));
DiagDecoh_qqp(:,iq,iqp) = PulseShape('gauss',[2*NDtP-1 FWHMqqp],1);
DiagDecoh_qqp(:,iq,iqp) = DiagDecoh_qqp(:,iq,iqp)./sum(DiagDecoh_qqp(:,iq,iqp));
end
end
rhoDtP = zeros(NDtP);
for iq = 1:NHH
for iqp = 1:NHH
rhot_qqp = sqrt(I_HH(iq))*sqrt(I_HH(iqp))*exp(-1i*2*pi*(E0_IR*HH(iq)-E0-Ip)*eV/h*nDtP(:)*dt -1i*(Phi(iq)-Phi(iqp))+1i*2*pi*(E0_IR*HH(iqp)-E0-Ip)*eV/h*nDtP*dt);
rhot_qqp = AddDiagDecoh(rhot_qqp,DiagDecoh_qqp(:,iq,iqp));
rhoDtP = rhoDtP + rhot_qqp;
end
end
rhoDtP = rhoDtP.*(shape(:)*shape);
rhoDtP = rhoDtP/trace(rhoDtP);
rhot = zeros(Nt_P);
rhot(Nzp/2+1:Nzp/2+NDtP,Nzp/2+1:Nzp/2+NDtP) = rhoDtP;
rhoE = SwitchStateRep(rhot,'Rhot2RhoE',1);
Purity(rhot)
sigma = eig(rhot);
figure(1414)
semilogy(abs(sigma),'ro-')
[V,D] = eigs(rhot,Nmodes);
P = V*sqrt(D);
RhoEig = P*P';
nu = Purity(RhoEig);
rhoEEig = SwitchStateRep(RhoEig,'Rhot2RhoE',1);
figure(1515)
imagesc((abs(rhoEEig)))
colormap(Carte_de_couleurs)
else
P = zeros(NDtP,1);
for iHH = 1:NHH
P = P + sqrt(I_HH(iHH))*exp(-1i*2*pi*(E0_IR*HH(iHH)-E0-Ip)*eV/h*nDtP(:)*dt);
end
P = P.*shape(:);
P = vertcat(zeros(Nzp/2,1),P,zeros(Nzp/2,1));
end
rhoE = SwitchStateRep(rhot,'Rhot2RhoE',1);
rhoDE=RhoE_2_RhoDE(rhoE,Ntau,Nc);
rhoE2=RhoDE_2_RhoE(rhoDE,Ntau,Nc);
%% Préallocation des matrices
E0_IR = 1.55;%Central photon energy of laser dressing pulse (eV)
T_IR = h/(E0_IR*eV);%Period of the dressing pulse (sec)
omega_IR = 2*pi/T_IR;
t = (-N/2:N/2-1)*dt;
p = KE2p(KE);%kinetic energy to momentum
Intensity = 8;%TW/cm2
IntensityW = Intensity*1e4*TW;%W/m2
FWHM_IR_intensity = 3*fs;
FWHM_IRfield = FWHM_IR_intensity*sqrt(2);
%Gaussian IR pulse
shape = PulseShape('gauss',[N round(FWHM_IRfield/dt)],1);
A_IR0 = -sin(omega_IR*t(:)+pi/2).*shape(:)/omega_IR;
A_IR = A_IR0*sqrt(2*IntensityW/(c*epsilon0));
%
% figure
% plot(t/fs,A_IR)
[Phi_mod,intA,intA2] = Phase_modulator(A_IR,dt,p);
Gt = exp(1i*Phi_mod);
%Crot =
% figure(1851)
% imagesc(t/fs,Ec,Phi_mod)
% title('Phase modulator \phi(\epsilon_c,t)')
% ylabel('Energy \epsilon_c (eV)')
% xlabel('time (fs)')
% colormap('jet')
Gtild = (ifft(Gt,[],2));
Gtild =fftshift(Gtild.',1);
Ef = Ec;
% figure(1852)
% imagesc(Ef,Ec,abs(Gtild))
% ylabel('\epsilon_1')
% xlabel('p')
% title('\tilde{G}(p,\epsilon_1)')
% colormap(Carte_de_couleurs)
%%
% 2) Pré-allocation CPU
Ctot = zeros(N, 2*NE, Ntau);
% 3) Boucle parallèle
parfor ipf = 1:N
% Produit externe sur CPU
g_vec = Gtild(:, ipf); % N×1
C = g_vec * g_vec.'; % N×N
% Passage en ρDE et rotation
Crot = fftshift( RhoE_2_RhoDE(C, Ntau, Nc, false), 2 ); % N×(NE?×Ntau?)
% Extension et décalage
Cext = [ Crot; zeros(NE, Ntau) ];
Csh = fftshift( circshift(Cext, [ipf-1, 0]) );
% Stockage slice-safe
Ctot(ipf, :, :) = Csh;
end
Stild = zeros(N,Ntau);
for iDE = 1:Ntau
vRhoDE = rhoDE(:,iDE);
% Quart haut : les premières 2 lignes
M1 =squeeze(Ctot(:,1:NE/2,iDE));
% Quart bas : les dernières 2 lignes
M2 = squeeze(Ctot(:,NE/2+NE+1:end,iDE));
size(M1)
size(M2)
% Mixer les quarts (alternance des lignes)
M =fftshift(horzcat(M1,M2),2);
% M =squeeze(Ctot(:,:,iDE));
% M = squeeze(Ctot(:,:,iDE));
% M_new=zeros(NE);
M_new=M;
% for i=1:2:2*N2
% for j=1:N2
% if mod(i - j, 2) == 0
% M_new(j,(i+1)/2) =M(j,i);
% else
% M_new(j,(i+1)/2) =M(j,i+1);
% end
% end
% end
% M_new=fftshift(M_new,1);
size(M_new)
size(vRhoDE)
Stild(:,iDE) = M_new*vRhoDE;
Stild(isnan(Stild)) = 0;
end
%% Ajout de bruit et préparation données
S = fft(Stild,[],2);
noise_level = 0e-6;
noise = (randn(size(S)) + 1i*randn(size(S)))/sqrt(2);
modS = abs(S);
phiS = angle(S);
modS_noisy = modS + noise_level * randn(size(S));
S_noisy = max(modS_noisy, 0) .* exp(1i * phiS);
Stild2 = ifft(S_noisy,[],2);
figure;
subplot(2,1,1);
imagesc(abs(S)); colorbar; title('Ideal S');
subplot(2,1,2);
imagesc(abs(Stild2)); colorbar; title('Noisy S');
figure;
subplot(2,1,1);
imagesc(angle(S_noisy)); colorbar; title('Ideal S');
subplot(2,1,2);
imagesc(angle(S)); colorbar; title('Noisy S');
%% Préparation des opérateurs pour la reconstruction
NE = 348; e = ones(NE,1);
D1 = spdiags([-e e],[0 1],NE,NE);
D2 = spdiags([e -2*e e],[-1 0 1],NE,NE);
W2 = D1'*D1; L2 = D2'*D2;
g1_list = logspace(-6,5,20); g2_list = logspace(-6,5,20); g3_list = logspace(-4,4,8);
iter_max=100
%=========================================================================
% Step 1: Tikh+TV Solve (to get the PHASE)
%=========================================================================
% Opérateur de gradient (TV 1D)
e = ones(NE,1);
D = spdiags([-e e], [0 1], NE, NE);
%--- 0) Load your data & set dimensions -------------------------------
% Ctot: NE × NE_total × Ntau
% Stild2: N × Ntau
% L2, W2: NE × NE
% NE, Ntau, Nc, N: scalars
% RhoDE_2_RhoE, ensurePositiveDefinite: function handles
%--- 1) Precompute finite-diff operator --------------------------------
e = ones(NE,1);
D1 = spdiags([-e e],[0 1],NE,NE);
DtD = D1' * D1;
var_d = var(Stild2, 0, 2); % N×1
%--- 2) Allocate outputs ------------------------------------------------
rhoDE_modulus = zeros(NE, Ntau);
phi_init = zeros(NE, Ntau);
gamma1_opts = zeros(Ntau,1);
gamma2_opts = zeros(Ntau,1);
%--- 3) First parfor: Tikhonov+TV solve to get module and raw phase -----
%====================================================================
% Parfor joint: Module by Tikh+TV + small-angle phase correction
%====================================================================

E = ones(NE,1);
D1 = spdiags([-E E], [0 1], NE-1, NE);
D2 = spdiags([E -2*E E], [-1 0 1], NE, NE);
alpha1 = 1e-4;
alpha2 = 1e-4;
Reg_block = alpha1 * blkdiag(D1'*D1, D1'*D1) + alpha2 * blkdiag(D2'*D2, D2'*D2);
%--- Boucle parfor ---------------------------------------------------
% Paramètres
g1_list = logspace(-4.5,5,20);
g2_list = logspace(-4.5,5,20);
alpha_tv = 0.0; % pas de TV pour le module
% Pré‐allocation
%--- Opérateurs de dérivées 1D -----------------------------------------
e = ones(NE,1);
D1 = spdiags([-e e],[0 1],NE,NE); % 1ʳᵉ dérivée (TV si besoin)
D2 = spdiags([ e -2*e e],[-1 0 1],NE,NE);% 2ᵉ dérivée (Laplacien 1D)
%--- Préallocation -----------------------------------------------------
rhoDE_modulus = zeros(NE, Ntau);
phi_init = zeros(NE, Ntau);
gamma1_opts = zeros(Ntau,1);
gamma2_opts = zeros(Ntau,1);
g3_list = logspace(-6, 4, 10); % 10 candidats pour la régularisation phase
lambda_lap_phase=1e-10;
P2 = D1' * D1; % pénalité sur variation de phase
gamma3_opts = zeros(Ntau,1);
parfor iDE = 1:Ntau
% a) Build M and original data
M1 = squeeze(Ctot(:,1:NE/2, iDE));
M2 = squeeze(Ctot(:,NE/2+NE+1:end, iDE));
M = fftshift([M1, M2], 2);
d = Stild2(:, iDE);
% b) Truncated‐SVD cleanup
[U,S2,V] = svd(M,'econ');
s = diag(S2);
r = max(find(s>=0.01*s(1),1,'last'),10);
Mclean = U(:,1:r)*S2(1:r,1:r)*V(:,1:r)';
% c) Linear solve for complex x
A_stack = [ real(Mclean), -imag(Mclean);
imag(Mclean), real(Mclean) ];
b_stack = [ real(d); imag(d) ];
x_sol = (A_stack'*A_stack + Reg_block) \ (A_stack'*b_stack);
x_cplx = x_sol(1:NE) + 1i*x_sol(NE+1:end);
% d) Median‑filter baseline removal on module
abs_m = abs(x_cplx);
win = max(5, round(NE/10)); % window ~NE/10
baseline = medfilt1(abs_m, win); % estimate slow background
m_ref = abs_m - baseline; % remove background
m_ref(m_ref < 0) = 0; % clamp negatives to zero
thr = 0.2 * max(m_ref); % 20% of the max
m_ref(m_ref < thr) = 0; % threshold weak residuals
% d) Store initial module and phase, plus sinogram prediction
rhoDE_modulus(:,iDE) = abs(m_ref);
phi_init(:,iDE) = angle(x_cplx);
end
% Phase finale
rhoDE_phase = phi_init;
%% Reconstruction finale
rhoDE_weighted = rhoDE_modulus .* exp(1i * phi_init);
rhoE_weighted = RhoDE_2_RhoE(rhoDE_weighted, Ntau, Nc);
% % 1) Hermitian symétrisation
% rhoE_rec = (rhoE_rec + rhoE_rec')/2;
% % 2) Décomposition spectrale
% [V,D] = eig(rhoE_rec);
% d = diag(D);
% % 3) Calcul du décalage minimal δ pour que toutes les λ_i+δ ≥ 0
% tol = 0;
% delta = max(0, -min(d) + tol);
% % 4) Reconstruction PSD en ne touchant qu'aux valeurs propres
% rhoE_pd = V * (D + delta*eye(size(D))) * V';
% % 5) Renormalisation de la trace
% rhoE_pd = rhoE_pd / trace(rhoE_pd);
% % 1) Calculez les valeurs propres de la matrice corrigée
% [V,D] = eig(rhoE_pd);
% eigvals = diag(D);
% % 2) Vérifiez qu’elles sont toutes ≥ 0 (à un petit tolérance près)
% tol_eig = -1e-12; % tolérance numérique
% if all(eigvals > tol_eig)
% disp('✅ Matrice semi‑définie positive (toutes les λ ≥ 0).');
% else
% warning('⚠️ Certaines valeurs propres sont négatives : min(λ) = %g', min(eigvals));
% end
% --- vous pouvez ensuite réaffecter, afficher, etc.
% rhoE_weighted = rhoE_pd;
figure(138);
subplot(1,2,1); imagesc(abs(rhoE_weighted)); axis image; colorbar; title('|ρ| final');
subplot(1,2,2); imagesc(angle(rhoE_weighted)); axis image; colorbar; title('∠ρ low-rank');
figure(139)
subplot(1,2,1); imagesc(abs(rhoE)); axis image; colorbar; title('|ρ| final');
subplot(1,2,2); imagesc(angle(rhoE)); axis image; colorbar; title('∠ρ low-rank');
%%
% Bootstrap par rééchantillonnage des résidus
Nboot = 20; % Nombre d’itérations bootstrap
% Initialisation stockage bootstrap
rhoDE_modulus_boot = zeros(NE, Ntau, Nboot);
phi_init_boot = zeros(NE, Ntau, Nboot);
% Matrices de dérivation
E = ones(NE,1);
D1 = spdiags([-E E], [0 1], NE-1, NE);
D2 = spdiags([E -2*E E], [-1 0 1], NE, NE);
alpha1 = 1e-4;
alpha2 = 1e-4;
Reg_block = alpha1 * blkdiag(D1'*D1, D1'*D1) + alpha2 * blkdiag(D2'*D2, D2'*D2);
% Reconstruction initiale complète (servira de base pour le bootstrap)
phi_init = zeros(NE, Ntau);
rhoDE_modulus = zeros(NE, Ntau);
d_hat_store = zeros(size(Stild2));
parfor iDE = 1:Ntau
% a) Build M and original data
M1 = squeeze(Ctot(:,1:NE/2, iDE));
M2 = squeeze(Ctot(:,NE/2+NE+1:end, iDE));
M = fftshift([M1, M2], 2);
d = Stild2(:, iDE);
% b) Truncated‐SVD cleanup
[U,S2,V] = svd(M,'econ');
s = diag(S2);
r = max(find(s>=0.01*s(1),1,'last'),10);
Mclean = U(:,1:r)*S2(1:r,1:r)*V(:,1:r)';
% c) Linear solve for complex x
A_stack = [ real(Mclean), -imag(Mclean);
imag(Mclean), real(Mclean) ];
b_stack = [ real(d); imag(d) ];
x_sol = (A_stack'*A_stack + Reg_block) \ (A_stack'*b_stack);
x_cplx = x_sol(1:NE) + 1i*x_sol(NE+1:end);
% d) Median‑filter baseline removal on module
abs_m = abs(x_cplx);
win = max(5, round(NE/10)); % window ~NE/10
baseline = medfilt1(abs_m, win); % estimate slow background
m_ref = abs_m - baseline; % remove background
m_ref(m_ref < 0) = 0; % clamp negatives to zero
thr = 0.2 * max(m_ref); % 20% of the max
m_ref(m_ref < thr) = 0;
% d) Store initial module and phase, plus sinogram prediction
rhoDE_modulus(:,iDE) = abs(x_cplx);
phi_init(:,iDE) = angle(x_cplx);
% d) stockage raw
rhoDE_raw(:,iDE) = x_cplx;
end
% % a) passage DE→E
% rhoE = RhoDE_2_RhoE(rhoDE_raw, Ntau, Nc);
% rhoE = (rhoE + rhoE')/2; % hermitian
%
% % b) décalage spectral minimal δ pour λ_min+δ≥0
% [Vb,Db] = eig(rhoE);
% db = diag(Db);
% sig1 = max(0, -min(db) + 1e-12);
% rhoE_pd = Vb * (Db + sig1*eye(size(Db))) * Vb';
% rhoE_pd = rhoE_pd / trace(rhoE_pd);
%
% % c) retour E→DE
% rhoDE_pd = RhoE_2_RhoDE(rhoE_pd, Ntau, Nc);
% 5) Maintenant vous pouvez construire d_hat_store ou tout autre post‑traitement
parfor iDE = 1:Ntau
M1 = squeeze(Ctot(:,1:NE/2, iDE));
M2 = squeeze(Ctot(:,NE/2+NE+1:end,iDE));
M = fftshift([M1, M2],2);
rhoDE_slice = rhoDE_raw(:,iDE);
d_hat_store(:,iDE) = M * rhoDE_slice;
end
% --- 2) Bootstrap sur résidus (module‐only + median‐baseline) ---
rhoDE_mod_boot = zeros(NE, Ntau, Nboot);
phi_boot = zeros(NE, Ntau, Nboot);
% Pré‑allocation
rhoDE_boot = zeros(NE, Ntau, Nboot); % sinogrammes DE bootstrap
rhoE_boot = zeros(250, 250, Nboot); % matrices densité EE bootstrap
for b = 1:Nboot
% --- 1) Reconstruction module+phase en DE pour chaque iDE ---
parfor iDE = 1:Ntau
M1 = squeeze(Ctot(:,1:NE/2, iDE));
M2 = squeeze(Ctot(:,NE/2+NE+1:end,iDE));
M = fftshift([M1, M2], 2);
d0 = Stild2(:, iDE);
d_hat = d_hat_store(:, iDE);
% residual bootstrap
res = d0 - d_hat;
d_star = d_hat + res(randsample(numel(res),numel(res),true));
% raw LS solve
[U,S2,V] = svd(M,'econ'); s=diag(S2);
r = max(find(s>=0.01*s(1),1,'last'),10);
Mclean = U(:,1:r)*S2(1:r,1:r)*V(:,1:r)';
A_stack = [ real(Mclean), -imag(Mclean);
imag(Mclean), real(Mclean) ];
b_stack = [ real(d_star); imag(d_star) ];
x_sol = (A_stack'*A_stack + Reg_block) \ (A_stack'*b_stack);
x_cplx = x_sol(1:NE) + 1i*x_sol(NE+1:end);
% median‑baseline removal on module
abs_m = abs(x_cplx);
win = max(5,round(NE/10));
baseline = medfilt1(abs_m,win);
m_ref = abs_m - baseline;
thr = 0.2*max(m_ref); m_ref(m_ref<thr)=0;
% keep original phase
phi_ref = angle(x_cplx);
% store sinogram DE
rhoDE_boot(:,iDE,b) = m_ref .* exp(1i*phi_ref);
end
% --- 2) Rotation DE→E + PSD projection, une seule fois par b ---
sinogram = rhoDE_boot(:,:,b); % NE×Ntau
rhoE_rec = RhoDE_2_RhoE(sinogram, Ntau, Nc); % Nc×Nc
% % Hermitian symétrisation
rhoE_rec = (rhoE_rec + rhoE_rec')/2;

% Spectral shift minimal pour PSD
[Vb,Db] = eig(rhoE_rec);
db = diag(Db);
tol = 1e-12;
delta = max(0, -min(db) + tol);
rhoE_pd = Vb * (Db + delta*eye(size(Db))) * Vb';
%
% % Renormalisation trace=1
rhoE_boot(:,:,b) = rhoE_pd;
[V,D] = eig(rhoE_boot(:,:,b));
eigvals = diag(D);
% 2) Vérifiez qu’elles sont toutes ≥ 0 (à un petit tolérance près)
tol_eig = -1e-12; % tolérance numérique
if all(eigvals > tol_eig)
disp('✅ Matrice semi‑définie positive (toutes les λ ≥ 0).');
else
warning('⚠️ Certaines valeurs propres sont négatives : min(λ) = %g', min(eigvals));
end
fprintf('Bootstrap %d/%d terminé\n', b, Nboot);
end
% --- 3) Post‑processing ---
rho_mod_mean = mean( rhoE_boot, 3 );
rho_mod_std = std( rhoE_boot, 0, 3 );
phi_mean = mean( phi_boot, 3 );
phi_std = std( phi_boot, 0, 3 );
rho_mod_var = var( rhoE_boot,0, 3 );
phi_var =var( phi_boot, 0, 3 );
% Moyenne et écart-type bootstrap
figure(130);
subplot(1,2,1); imagesc(abs(rho_mod_mean)); axis image; colorbar; title('|ρ| final');
subplot(1,2,2); imagesc(angle(rho_mod_mean)); axis image; colorbar; title('∠ρ low-rank');
figure(131);
subplot(1,2,1); imagesc(mod_var); axis image; colorbar; title('|ρ| final');
subplot(1,2,2); imagesc(abs(rhoE2)); axis image; colorbar; title('∠ρ low-rank');