function rhoE = RhoDE_2_RhoE(rhoDE,Ntau,Nc)

NE = Sampling_Sami(Nc,Ntau);

Ni = floor((Ntau-1)/2);


rhoErot_full = zeros(2*NE-1);
rhoErot_full(Ni+1:Ni+2*Nc,NE-(Ntau-1)/2:NE+(Ntau-1)/2) = rhoDE;



rhoErot_fullfull = zeros(2*NE-1);
for ic = 1:2*NE-1
    rhoErot_fullfull(:,ic) = circshift(rhoErot_full(:,ic),-ic+NE);
end


figure(819)
imagesc(abs(rhoErot_fullfull))%,[0 max(vdiag)])
axis equal
colormap(Carte_de_couleurs)

rhoE = zeros(NE);
for ir = 1:NE
    V = NE-ir+1 + (0:NE-1);
    rhoE(ir,:) = rhoErot_fullfull(2*ir-1,V);
end

figure(8181)
subplot(1,3,1)
imagesc(abs(rhoErot_full))
axis equal
%colormap(Carte_de_couleurs)

subplot(1,3,2)
imagesc(abs(rhoErot_fullfull))%,[0 max(vdiag)])
axis equal
%colormap(Carte_de_couleurs)


subplot(1,3,3)
imagesc(abs(rhoE))
axis equal
%colormap(Carte_de_couleurs)
