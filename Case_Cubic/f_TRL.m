% Reference paper: Ying, Li & Fu, Generalized synthetic inflow generation %
% method for divergence-free inhomogeneous turbulence, JCP, 2025          %
% Authors: Anjia Ying, Zhigang Li, Lin Fu                                 %

%%%%% Function: Construct the 3D covariance matrix of homogeneous 
%%%%% turbulence with the TRL model

function data = f_TRL(z,t,U,nu,Ruvw,mL)

    d_z = z(2)-z(1);
    d_t = t(2)-t(1);
    
    Nz = length(z);
    Nt = length(t);
    
    UU = Ruvw(1,1);
    VV = Ruvw(2,2);
    WW = Ruvw(3,3);
    
    L1 = mL(1);
    L2 = mL(2);
    L3 = mL(3);
    
    %% Wavenumbers and frequencies corresponding to the grid and temporal settings
    Nkz = Nz;
    Nft = Nt;
    
    % Nyquist frequency
    Nq_z = 2*pi/d_z/2;
    Nq_t = 2*pi/d_t/2;
    kz0= 2*Nq_z*(0:Nkz-1)'/Nkz;
    f0 = 2*Nq_t*(0:Nft-1)'/Nft;
    
    % Frequencies from 0 to 2Nq
    index_z = find(kz0 >= Nq_z*0.9999,1)-1;
    index_t = find(f0  >= Nq_t*0.9999,1)-1;
    
    kz = kz0;
    kz(index_z+1:Nkz) = kz0(index_z+1:Nkz)-2*Nq_z;
    dkz = kz(2)-kz(1);
    
    ft = f0;
    ft(index_t+1:Nft)  = f0(index_t+1:Nft) -2*Nq_t;
    dft = ft(2)-ft(1);
    
    Nz_mod = Nz/2;
    Nt_mod = Nft/2;
    
    % Equilibrim streamwise wavenumber kxm (when determining the spectra, the values of wavenumbers are set to be positive)
    ftm = max(ft(1:Nt_mod),1e-7);
    Nftm = length(ftm);
    kxm  = ftm/U;
    
    %% E_iso (1D) from Pope's spectra model
    K = (UU+VV+WW) /2;
    E_iso = f_Spectra_Pope(K,L1',nu , kxm);  % Nkx*Nkx*Nkx * Ny
    
    %% Interpolate the E_HAT (3D) from E_iso (1D)
    sigma_L2 = L2/L1;
    sigma_L3 = L3/L1;
    
    kk3 = max(sqrt(kxm.^2 + (kz.*sigma_L2)'.^2 + permute((kz.*sigma_L3),[3 2 1]).^2),1e-6);
    kk3_1d = kk3(:);
    E_HAT = reshape( interp1(kxm , E_iso./(4*pi*kxm.^2) , kk3_1d) , [Nftm,Nkz,Nkz] );   % Nkx,Nkz,Nkz
    E_HAT(1,:,:) = 0.5*E_HAT(1,:,:);
    
    E_HAT(find(isnan(E_HAT))) = 0;
    E_HAT = E_HAT ./ sum(E_HAT*dft*dkz^2,[1 2 3]);  % Ensure the unit energy for UVW
    
    %% Scale E_HAT to the null space of wavenumbers
    scale_C = zeros(Nftm,Nkz,Nkz,3,3);
    k1 = kxm;
    k2 = kz';
    k3 = permute( kz ,[3 2 1]);
    kkn = max(sqrt(k1.^2 + k2.^2 + k3.^2),1e-6);
    
    scale_11 = 1 - ( k1.* k1 )./kkn.^2 ;
    scale_22 = 1 - ( k2.* k2 )./kkn.^2 ;
    scale_33 = 1 - ( k3.* k3 )./kkn.^2 ;
    
    scale_12 = -( k1.* k2 )./kkn.^2 ;
    scale_13 = -( k1.* k3 )./kkn.^2 ;
    scale_23 = -( k2.* k3 )./kkn.^2 ;
    
    scale_C(:,:,:,1,1) = scale_11;
    scale_C(:,:,:,2,2) = scale_22;
    scale_C(:,:,:,3,3) = scale_33;
    
    scale_C(:,:,:,1,2) = scale_12;
    scale_C(:,:,:,1,3) = scale_13;
    scale_C(:,:,:,2,3) = scale_23;
    
    scale_C(:,:,:,2,1) = scale_12;
    scale_C(:,:,:,3,1) = scale_13;
    scale_C(:,:,:,3,2) = scale_23;
    
    Cs_HAT = real(scale_C .* E_HAT); % Nkx*Nkx*Nkx * 3 * 3
    
    %% Normalize the spectrum
    uu = squeeze(1./mean(Cs_HAT(:,:,:,1,1),[1 2 3]));
    vv = squeeze(1./mean(Cs_HAT(:,:,:,2,2),[1 2 3]));
    ww = squeeze(1./mean(Cs_HAT(:,:,:,3,3),[1 2 3]));
    uvw_sc = sqrt([uu,vv,ww]);
    uvw_sc3 = permute(uvw_sc,[5 4 3 2 1]) .* permute(uvw_sc,[5 4 3 1 2]);
    Cs_HAT_norm = Cs_HAT .* uvw_sc3;
    
    %% Match the target Reynolds stress
    RL = chol(Ruvw);
    Rij = RL;
    
    CST_H3 = permute(Cs_HAT_norm,[4 5 1 2 3]);
    CST_temp1 = pagemtimes(Rij',CST_H3);
    CST_temp2 = pagemtimes(CST_temp1,Rij);
    
    % Rescale the Reynods stress value
    Euu_rec = sum(CST_temp2(1,1,:,:)*dft*dkz,[3 4]);
    
    factorR = max(UU)/max(Euu_rec);
    
    CST_temp3 = real(factorR*CST_temp2);
    
    CST = 4*CST_temp3(:,:,:,1:Nz_mod,1:Nz_mod);
    CST(1,3,:,:,:) = 0;
    CST(3,1,:,:,:) = 0;
    CST(2,3,:,:,:) = 0;
    CST(3,2,:,:,:) = 0;
    CST(:,:,:,1,:) = 0.5*CST(:,:,:,1,:);
    CST(:,:,:,:,1) = 0.5*CST(:,:,:,:,1);
    
    data.ft = ft;
    data.kz = kz;
    data.Nt_mod = Nt_mod;
    data.Nz_mod = Nz_mod;
    data.CST = CST;

end
