function [PsiChan1,PsiChan2,rhoc] = DensMatChanSplit(Rho,indxSplit,Niter,display)
%Function that implements the two-channel separation of a density matrix
%composed of overlapping channels
%
% USE: [PsiChan1,PsiChan2,Err] = DensMatChanSplit(Rho,indxSplit,Niter)
%                                       OR
%      [PsiChan1,PsiChan2,Err] = DensMatChanSplit(Rho,indxSplit,Niter,display)
%
% INPUTS:
% Rho: complex matrix of size (NxN) representing the density matrix with
%      overlapping channels
%
% indxSplit: index of the energy below which the 2nd channel is zero. It can be :
%            - an integer between 1 and N, then the state vector of the 2nd
%            channel is constrained to be zero below indxSplit
%            - a two-element vector [min(indxSplit) max(indxSplit)] containing  the range of indices to
%            scan in order to find the best threshold index
%
% Niter: integer corresponding to the number of iteration of the algorithm
%
% display (optional): logical (0 or 1), turns on/off display
%
% OUTPUTS:
% PsiChan1,PsiChan2: complex vectors (Nx1) containing the vectors of both channels
%
% Err : Evolution of the Error (NRMSE) between the input and the fitted
%       density matrices as a function of the iteration number. It can be:
%       - a vector of size (1xNiter) if indxSplit is an integer
%       - a matrix of size (NindxSplit x Niter) if indxSplit is a
%       two-element vector, in that case NindxSplit is the number of scanned values, i.e. the size of min(indxSplit):max(indxSplit)
%
% C. Bourassin-Bouchet, 13/01/2025


if nargin == 3
    display = 0;
end

if numel(indxSplit) == 1
    splittype = 'manual';
    
    indxSplitmin = indxSplit;
    indxSplitmax = indxSplit;
    
elseif numel(indxSplit) == 2
    splittype = 'auto';
    
    indxSplitmin = min(indxSplit);
    indxSplitmax = max(indxSplit);
end

N = size(Rho,1);
NindxSplit = length(indxSplitmin:indxSplitmax);

Err = zeros(NindxSplit,Niter);


%% Splitting algorithm
cmpt_Mask = 0;
for iMask = indxSplitmin:indxSplitmax
    cmpt_Mask = cmpt_Mask+1;
    
    %Random initial guess
    PSI1 = rand(N,5)+1i*rand(N,5);
    
    Wrtvd = ones(N,5);
    Wrtvd(1:iMask,2) = 0;%Channel 2 is set to zero below the iMask index
    Wrtvd(1:iMask,3) = 0;
    Wrtvd(1:iMask,4) = 0;
    Wrtvd(1:iMask,5) = 0;
    %Weighted low-rank approximation algorithm (weighted subspace iteration)
    for iter = 1:Niter
        PSI2 = PSI1\Rho;
        
        PSI2 = (PSI2').*Wrtvd;
        
        PSI1 = PSI2\Rho;
        
        PSI1 = (PSI1').*Wrtvd;
        
        %NRMSE between the original and fitted density matrix 
        Err(cmpt_Mask,iter) = NRMSError(Rho,PSI1*PSI2');
    end
end


%% Finding the best threshold index
if strcmp(splittype,'auto')
    
    [~,OptimSplitIndx] = min(Err(:,Niter));

    PSI1 = rand(N,5)+1i*rand(N,5);
    Wrtvd = ones(N,5);
    Wrtvd(1:(indxSplitmin+OptimSplitIndx-1),2) = 0;
    Wrtvd(1:(indxSplitmin+OptimSplitIndx-1),3) = 0;
    Wrtvd(1:(indxSplitmin+OptimSplitIndx-1),4) = 0;
    Wrtvd(1:(indxSplitmin+OptimSplitIndx-1),5) = 0;
    
    for iter = 1:Niter
        PSI2 = PSI1\Rho;
        
        PSI2 = (PSI2').*Wrtvd;
        
        PSI1 = PSI2\Rho;
        
        PSI1 = (PSI1').*Wrtvd;
    end
end


%% Final vector estimation

RhoChan1 = PSI1(:,1)*PSI2(:,1)';
RhoChan2 = PSI1(:,2)*PSI2(:,2)';

[V,D] = eigs(RhoChan1,1);
PsiChan1 = sqrt(D)*V;

[V,D] = eigs(RhoChan2,1);
PsiChan2 = sqrt(D)*V;

rhoc = PSI1*PSI2';

%% Optional display
if display
    
    n = 0:N-1;
    
    figure
    subplot(2,4,1)
    imagesc(abs(Rho))
    xlabel('Energy (index)')
    ylabel('Energy^prime (index)')
    title('Original |\rho|')
    colormap(Carte_de_couleurs)
    
    subplot(2,4,2)
    imagesc(angle(Rho))
    xlabel('Energy (index)')
    ylabel('Energy^prime (index)')
    title('Original arg(\rho)')
    colormap(Carte_de_couleurs)
    
    subplot(2,4,3)
    plot(n,abs(PsiChan1),'ko-',n,angle(PsiChan1),'ro-')
    xlabel('Energy (index)')
    legend('Modulus','Phase')
    title('Retrieved \Psi_1')
    
    
    if strcmp(splittype,'auto')

        subplot(2,4,4)
        semilogy(indxSplitmin:indxSplitmax,Err(:,Niter),'ro-')
        grid on
        ylabel('NRMSE')
        xlabel('Threshold index')
        title('Final NRSME')
        
        subplot(2,4,8)
        semilogy(1:Niter,Err(OptimSplitIndx,:),'ro-')
        grid on
        ylabel('NRMSE')
        xlabel('iteration number')
        title('NRSME @optimal threshold index')
        
    else
        
        subplot(2,4,[4 8])
        semilogy(1:Niter,Err,'ro-')
        grid on
        ylabel('NRMSE')
        xlabel('iteration number')
        
    end
    
    
    subplot(2,4,5)
    imagesc(abs(PSI1*PSI2'))
    xlabel('Energy (index)')
    ylabel('Energy^prime (index)')
    title('Fitted |\rho|')
    colormap(Carte_de_couleurs)
    
    subplot(2,4,6)
    imagesc(angle(PSI1*PSI2'))
    xlabel('Energy (index)')
    ylabel('Energy^prime (index)')
    title('Fitted arg(\rho)')
    colormap(Carte_de_couleurs)
    
    subplot(2,4,7)
    plot(n,abs(PsiChan2),'ko-',n,angle(PsiChan2),'ro-',n,Wrtvd(:,2),'bo-')
    xlabel('Energy (index)')
    legend('Modulus','Phase','Threshold')
    title('Retrieved \Psi_2')
    
end

