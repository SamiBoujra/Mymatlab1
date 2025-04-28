
close all
clc
clear



%% Constants
load constants.mat



%% Ef sampling

Nf = 256;
nf = 0:Nf-1;
nfc= (-Nf/2:Nf/2-1);

E0f = 10;%eV
DeltaEf = 80;
dEf = DeltaEf/Nf;

dt = h/(DeltaEf*eV);%time step (sec)

t = (-Nf/2:Nf/2-1)*dt;

Ef = E0f + (0:Nf-1)*dEf;


%% delay/Delta E sampling

Ntau = 400;
dtau = 0.01*fs;
ntau = 0:Ntau-1;
tau = ntau*dtau;

dDeltaE = h/(Ntau*dtau*eV);
DeltaE = (-Ntau/2:Ntau/2-1)*dDeltaE;


%% Ec sampling

Z = 1;
Nc = round(2*Z*DeltaEf*eV/h*Ntau*dtau);



dEc = DeltaEf/Nc;
E0c = E0f;

nc = 0:Nc-1;
Ec = E0c + nc*dEc;



%% Density matrix

shape = PulseShape('sin2',[Nc Nc]);


NE = Sampling_Sami(Nc,Ntau);

%rhoE = shape(:)*shape;

ndiag = -NE:NE;
vdiag = abs(ndiag);
vdiag(vdiag == 0) = 6;
B = ones(NE,1)*vdiag;
M = zeros(NE);
rhoE = full(spdiags(B,ndiag,M));



figure(815)
imagesc(rhoE,[0 max(vdiag)])
colormap(Carte_de_couleurs)

%%
rhoDE = RhoE_2_RhoDE(rhoE,Ntau,Nc)


%%
rhoE2 = RhoDE_2_RhoE(rhoDE,Ntau,Nc)

%%
%% KE sampling

N = 256;
n = 0:N-1;
nc= (-N/2:N/2-1);

E0 = 10;%eV
DeltaE = 80;
dE = DeltaE/N;

dt = h/(DeltaE*eV);%time step (sec)
dE = DeltaE/N;
%% delay/Delta E sampling

Ntau = 150;
dtau = 0.01.*fs;
ntau = 0:Ntau-1;
tau = ntau*dtau;
dDeltaE = h_bar/(Ntau*dtau*eV);
DeltaE = (-Ntau/2:Ntau/2-1)*dDeltaE;
% Initialise Z et Nc
zCandidate = 0;
NcCandidate = 0;

% Tant que la condition est vérifiée, on incrémente zCandidate
continueLoop = true;

while continueLoop
    % On teste le prochain Z
    testZ = zCandidate + 1;
    testNc = round(2 * testZ * DeltaEf * (eV/h) * Ntau * dtau);

    if testNc <= N
        % Ce Z est valide, on met à jour
        zCandidate = testZ;
        NcCandidate = testNc;
    else
        % On a dépassé la limite, donc on arrête la boucle
        continueLoop = false;
    end
end



Z = max(0, zCandidate - 1);        % on évite Z < 0
Nc = round(2 * Z * DeltaEf * (eV/h) * Ntau * dtau);
NE =2*Nc;
t = (-N/2:N/2-1)*dt;

KE = E0 + (0:N-1)*dE;

shape = PulseShape('sin2',[N N/60]);

%First wavepacket
E0_EWP1 = 40;%central energy (eV)
P = shape(:).*exp(-1i*2*pi*(E0_EWP1-E0)*eV/h*n(:)*dt);

%Second wavepacket
E0_EWP2 = 20;%central energy (eV)
P(:,2) = 0.5*shape(:).*exp(-1i*2*pi*(E0_EWP2-E0)*eV/h*n(:)*dt);


[rhot,nu] = Density_matrix(P,1);
rhot = fftshift(rhot);
rhoE = SwitchStateRep(rhot,'Rhot2RhoE',1);
rhose=rhoE*rhoE'
figure(7)
imagesc(abs(rhoE2))
% Supposons que vous vouliez ajouter p zéros de chaque côté :
p = 100;

% Zero-padding sur les 4 côtés :
rhoE_padded = padarray(rhoE, [p p], 1, 'both');
rhoDE = RhoE_2_RhoDE(rhoE,Ntau,Nc);
rhoE2 = RhoDE_2_RhoE(rhoDE,Ntau,Nc);
Niter = 900;
indxSplit = 37;
[PsiChan1,PsiChan2,rhoc] = DensMatChanSplit(rhoE2,[25 50],Niter,1);
Rhos=PsiChan1*PsiChan2'; 
% rhoE2_bin(rhoE2 ~= 0) = 1;
figure(78)
imagesc(abs(rhoE))
figure(72)
imagesc(abs(Rhos))
figure(75)
imagesc(abs(rhoc))

%%
% A0 = 3;
% E0_IR = 1.55;%Central photon energy of laser dressing pulse (eV)
% T_IR = h/(E0_IR*eV);%Period of the dressing pulse (sec)
% dtG=0.08;
% df = 1/(N*dt);
% fIR = 1/T_IR;
% kfIR = 4;
% E0_IR = 1.55;%Central photon energy of laser dressing pulse (eV)
% T_IR = h/(E0_IR*eV*fs);%Period of the dressing pulse (sec)
% 
% sigt = 3;
% Phi_mod = A0.*exp(-(nc*dtG).^2/(2*sigt^2)).*cos(2*pi*nc*dtG/T_IR);
% Gt = fftshift(exp(1i*Phi_mod));

% figure(9184)
% plot(Phi_mod)

E0_IR = 1.55;%Central photon energy of laser dressing pulse (eV)
T_IR = h/(E0_IR*eV);%Period of the dressing pulse (sec)
omega_IR = 2*pi/T_IR;



p = KE2p(KE);%kinetic energy to momentum

Intensity = 8;%TW/cm2
IntensityW = Intensity*1e4*TW;%W/m2


FWHM_IR_intensity = 3*fs;
FWHM_IRfield = FWHM_IR_intensity*sqrt(2);

%Gaussian IR pulse
shape = PulseShape('gauss',[N round(FWHM_IRfield/dt)],1);

A_IR0 = -sin(omega_IR*t(:)+pi/2).*shape(:)/omega_IR;

A_IR = A_IR0*sqrt(2*IntensityW/(c*epsilon0));

figure
plot(t/fs,A_IR)

[Phi_mod,intA,intA2] = Phase_modulator(A_IR,dt,p);
Gt = exp(1i*Phi_mod);

% ish = 500;
% %Gtilde = circshift(ifft(Gt).',ish);
% Gtilde = 1/N*fft(Gt(:));
% C = fftshift(Gtilde*Gtilde');

%Crot = 
figure(1851)
imagesc(t/fs,Ec,Phi_mod)
title('Phase modulator \phi(\epsilon_c,t)')
ylabel('Energy \epsilon_c (eV)')
xlabel('time (fs)')
colormap('jet')

Gtild = (ifft(Gt,[],2));
Gtild =fftshift(Gtild.',1);
Ef = Ec;

figure(1852)
imagesc(Ef,Ec,abs(Gtild))
ylabel('\epsilon_1')
xlabel('p')
title('\tilde{G}(p,\epsilon_1)')
colormap(Carte_de_couleurs)


% figure(9185)
% subplot(2,1,1)
% plot(Ef,(abs(Gtilde)))
% subplot(2,1,2)
% imagesc(Ef,Ef,(abs(C)))
% colormap(Carte_de_couleurs)

%%

% Crot = RhoE_2_RhoDE(C,Ntau,Nc);
% N2=size(Crot,1);
% figure(1414)
% imagesc(abs(Crot.'))
% colormap(Carte_de_couleurs)

%%

Ctot = zeros(N,2*NE,Ntau);


for ipf = 1:N
    
    C = zeros(2*NE);
    C(1:N,1:N) = (Gtild(:,ipf)*Gtild(:,ipf)');
    
    Crot =fftshift(RhoE_2_RhoDE(C,Ntau,Nc,false),2);
    Crot2 = circshift(vertcat(Crot,zeros(NE,Ntau)),[ipf-1 0]);
    Crot2= fftshift(Crot2);
    
    Ctot(ipf,:,:) = Crot2(:,:);

     % figure(765)
     % imagesc(squeeze(abs(Ctot(ipf,:,:))))
     % colormap(Carte_de_couleurs)
     % pause(0.1)
end
%%

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
    %  for j=1:N2
    %     if mod(i - j, 2) == 0
    %         M_new(j,(i+1)/2) =M(j,i);
    %     else
    %         M_new(j,(i+1)/2) =M(j,i+1);
    %     end
    %  end
    % end
    % M_new=fftshift(M_new,1);
    size(M_new)
    size(vRhoDE)
    Stild(:,iDE) = M_new*vRhoDE;
% %     
      RhoE_rot_disp = rhoDE;
      RhoE_rot_disp(:,iDE) = max(max(rhoDE));
    
    figure(181)
    subplot(1,2,1)
    imagesc(abs(M_new))
    ylabel('E_f (pix)')
    xlabel('E_c (pix)')
    colormap(Carte_de_couleurs)
    
    subplot(1,2,2)
    imagesc(abs(RhoE_rot_disp))
    colormap(Carte_de_couleurs)
    pause(0.1)
    
end

Stild(isnan(Stild)) = 0;

%Stild(:,N/2+2:N) = conj(fliplr(Stild(:,2:N/2)));

figure(493)
plot((sum(abs(Stild))),'ro-')

figure(494)
semilogy((sum(abs(Stild))),'ro-')

figure(195)
imagesc(abs(Stild))
colormap(Carte_de_couleurs)
%%
S = fft(Stild,[],2);

figure(1958)
imagesc(abs(S))
colorbar
colormap(Carte_de_couleurs)
figure(1959)
plot(real(S(181,:)),'ro-')
%%
    % 2) Add noise to S (S_ideal).
    %    For example: Gaussian noise with SNR ~ 20 dB
    % -----------------------------------------------------------
 mean(abs(S));
 noise_level = 0*1e-6; 
 S_noisy = S + noise_level * randn(size(S));

 S_noisy1=ifft(S_noisy,[],2);
    % Optional: Poisson-based noise
    % S_noisy = poissrnd(max(0,S_ideal));

figure; 
subplot(2,1,1);
imagesc(abs(S)); colorbar; title('Ideal S');
subplot(2,1,2);
imagesc(abs(S_noisy)); colorbar; title('Noisy S');
Stild2 = ifft(S_noisy, [], 2);
% Preallocate
[N, twoNE, Ntau] = size(Ctot);
% twoNE should be 2*NE
NE_half = (twoNE / 2);  % = NE

% ----------------------------
% Reconstruction par GCV slice‑by‑slice
% ----------------------------
%% Pré‑requis : Ec, tau, Ctot, Stild2, E0_EWP1, E0_EWP2, dE, fs, N*1e6(

NE    = size(Ctot,2)/2;    % Ctot est N x (2*NE) x Ntau
Ntau  = size(Ctot,3);

% Paramètres pour le balayage gamma
gmin        = 1e-8;
gmax        = 1e-2;
ngamma      = 50;
gamma_list  = logspace(log10(gmin), log10(gmax), ngamma);
% build a diagonal weight matrix W for energy‑weighted Tikhonov
% e.g., penalize far‑off energies more heavily
KE1 = E0 + (0:NE-1)*dE;
E0_EPW1A=floor(E0_EWP1/0.73);
E0_EWP2A=floor(E0_EWP2/0.73);
% Paramètres de la gaussienne
% define your ROI edges
Emax =E0_EPW1A;
Emin = E0_EWP2A;

% build a logical mask: 1 inside [Emin,Emax], 0 outside
maskE = ( KE1 >= Emin-5 ) & ( KE1 <= Emax+5 );   
maskE=maskE.';% length NE

% if you want to allow two disjoint windows you can do, e.g.
% maskE = ( abs(KE1-E1)<=Δ1 ) | ( abs(KE1-E2)<=Δ2 );



% now form B = M'*M, Bd = M'*d as before...
% and solve    (B + γ·W)\Bd
figure(78000)
imagesc(tau, KE1, maskE);  % Affichage avec tau (temps) sur l'axe X et KE1 (énergie) sur l'axe Y
colorbar;
xlabel('Delay (fs)');
ylabel('Energy (eV)');
title('Template 2D-Gaussian Circular');

%% Reconstruction avec weighted Tikhonov
rhoDE_weighted = zeros(NE, Ntau);
for iDE=1:Ntau
% Boucle sur chaque tranche de délaior iDE = 1:Ntau
    % Construction de la matrice modèle M et vecteur data d
    M1 = squeeze(Ctot(:,1:NE/2,        iDE));
    M2 = squeeze(Ctot(:,NE/2+NE+1:end, iDE));
    M  = fftshift([M1, M2], 2);  % N×NE
    d  = Stild2(:,iDE);          % N×1

    
    w = double(~maskE);    % w(i)=0 for i in ROI, 1 outside
    W = diag( w );         % NE×NE penalty matrix
    % Pré-calculs pour Tikhonov
    B  = M' * M;                 % NE×NE
    Bd = M' * d;                 % NE×1

    % L-curve: calcul des normes pour chaque gamma
    Rnorm = zeros(ngamma,1);
    Xnorm = zeros(ngamma,1);
    for k = 1:ngamma
        A_k = B + gamma_list(k) * W;
        xk  = A_k \ Bd;
        Rnorm(k) = norm(M*xk - d);
        Xnorm(k) = norm(xk);
    end

    % Sélection automatique de gamma au "coin" de la L-curve
    gamma_opt = pick_lcurve_corner(Rnorm, Xnorm, gamma_list);

    % Résolution finale avec gamma optimal
    A_opt      = B + gamma_opt * W;
    x_final    = A_opt \ Bd;
    rhoDE_weighted(:, iDE) = x_final;

    fprintf('Slice %d: gamma_opt = %.2e\n', iDE, gamma_opt);
end

% Conversion en rhoE (fonction utilisateur)
rhoE_weighted = RhoDE_2_RhoE(rhoDE_weighted, Ntau, Nc);

% Affichage des résultats
figure(2);
subplot(1,2,1);
imagesc(tau, KE1, W);
title('Template 2D-Gaussian ROI');
xlabel('Delay (fs)'); ylabel('Energy (eV)');
colormap('jet');

subplot(1,2,2);
imagesc(tau, Ec, abs(rhoE_weighted));
title('Reconstruction weighted Tikhonov');
xlabel('Delay (fs)'); ylabel('Energy (eV)');
colormap('jet');

%% Fonctions auxiliaires
function gamma_opt = pick_lcurve_corner(Rnorm, Xnorm, gamma_list)
% pick_lcurve_corner  Trouve gamma au coin de la L-curve via courbure discrète

logR = log(Rnorm(:));
logX = log(Xnorm(:));
n = numel(logR);

% Tri par logR croissant pour définir un chemin
[~, order] = sort(logR);
Lr = logR(order);
Lx = logX(order);

% Dérivées par différences finies
d1x = diff(Lx) ./ diff(Lr);
d2x = diff(d1x) ./ diff(Lr(2:end));

% Calcul de la courbure κ_k = |x' y'' - y' x''| / ( (x'^2 + y'^2)^(3/2) )
kappa = nan(n,1);
for k = 2:n-1
    dx  = d1x(k-1);
    dr  = 1;  % dLr/dLr = 1
    ddx = d2x(k-1);
    % d2r = 0 par structure
    num = abs(dx*0 - dr*ddx);
    den = (dx^2 + dr^2)^(3/2);
    kappa(order(k)) = num/den;
end

% Choix de l'indice de courbure maximale
[~, idx] = max(kappa);
% Choix de l'indice de courbure maximale\,[~, idx] = max(kappa);
gamma_opt = gamma_list(idx);
end

figure
imagesc(abs(rhoDE_weighted))
%% 
%% 


%%
% Now visualize:
% rhoE3 = RhoDE_2_RhoE(rhoE_ls,Ntau,Nc);
W = ones(size(rhoE));
W=RhoE_2_RhoDE(W,Ntau,Nc);
W=RhoDE_2_RhoE(W,Ntau,Nc);

[PsiChan1,PsiChan2,rhoc] = DensMatChanSplit(rhoE_weighted ,[25 50],100);

rhoc_sm = smooth2D_with_smooth(rhoc, 40, 'moving');
figure; 
subplot(1,2,1); imagesc(abs(rhoc)); colorbar; title('Global');
subplot(1,2,2); imagesc(abs(rhoE2)); colorbar; title('Local');
    % 2) Add noise to S (S_ideal).
    %    For example: Gaussian noise with SNR ~ 20 dB
    % -----------------------------------------------------------
%  noise_level = 0.5*1e-6; 
%  S_noisy = S + noise_level * randn(size(S));
% S_noisy1=ifft(S_noisy,[],2)
%     % Optional: Poisson-based noise
%     % S_noisy = poissrnd(max(0,S_ideal));
% 
% figure; 
% subplot(2,1,1);
% imagesc(abs(S)); colorbar; title('Ideal S');
% subplot(2,1,2);
% imagesc(abs(S_noisy)); colorbar; title('Noisy S');
% Stild = ifft(S_noisy, [], 2);
% % Preallocate
% [N, twoNE, Ntau] = size(Ctot);
% % twoNE should be 2*NE
% NE_half = (twoNE / 2);  % = NE
% M_new   = zeros(N, NE, Ntau);
% 
% for iDE = 1:Ntau
%     M1 = squeeze(Ctot(:, 1:NE_half/2,         iDE));  % top half
%     M2 = squeeze(Ctot(:, NE_half/2+NE_half+1:end, iDE));  % bottom half
%     M  = [M1, M2];              % horizontally concatenate => (N x NE)
%     M  = fftshift(M, 2);        % shift columns
%     M_new(:,:,iDE) = M;
% end
% % Suppose we have:
% %   M_new{iDE}  - cell array of forward matrices
% %   S_noisy     - measured data (size [N x Ntau])
% %   rhoDE_true  - optional ground truth, same [NE x Ntau]
% [rhoGlob, rhoLoc] = two_step_wiener_min(M_new, S_noisy1);
% % Now visualize:
% rhoE3 = RhoDE_2_RhoE(rhoGlob,Ntau,Nc);
% [PsiChan1,PsiChan2,rhoc]=DensMatChanSplitRect(rhoGlob,Ntau,50,true);
% 
% figure; 
% subplot(1,2,1); imagesc(abs(rhoc)); colorbar; title('Global');
% subplot(1,2,2); imagesc(abs(rhoE)); colorbar; title('Local');