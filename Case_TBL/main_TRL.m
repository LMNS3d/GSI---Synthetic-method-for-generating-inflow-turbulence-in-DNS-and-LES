% Reference paper: Ying, Li & Fu, Generalized synthetic inflow generation %
% method for divergence-free inhomogeneous turbulence, JCP, 2025          %
% Authors: Anjia Ying, Zhigang Li, Lin Fu                                 %

%%%% Example program: Construct the 3D covariance matrix of TBL with the TRL model

clear;close all;

% Grid settings (consistent with the CFD computation)

Nt = 400;
Ny = 151;
Nz = 150;

Lx = 200;
Lz = 15;

nu = 1/1000;    % Kinetic viscosity (Be consistent with CFD computation)
ratio_X = 1;

% Mean quantities of velocities
load('./Mean_U_V.mat','Umean','Vmean','UU','VV','WW','UV')

Uinf = Umean(end);

t = (1:Nt)/Nt*Lx / Uinf;
z = (1:Nz)/Nz*Lz;

y = load('./yp.dat');

d_y = (y(3:end)-y(1:end-2))/2;
d_y = [d_y(1);d_y;d_y(end)];

d_t = t(2)-t(1);
d_z = z(2)-z(1);

Nkz = Nz;                    % z
Nft = Nt;

Nq_z = 2*pi/d_z/2;    % z
Nq_t = 2*pi/d_t/2;
kz0= 2*Nq_z*(0:Nkz-1)'/Nkz;
f0 = 2*Nq_t*(0:Nft-1)'/Nft;

index_z = find(kz0 > Nq_z,1)-1;
index_t = find(f0  > Nq_t,1)-1;

kz = kz0;
kz(index_z+1:Nkz) = kz0(index_z+1:Nkz)-2*Nq_z;

ft = f0;
ft(index_t+1:Nft)  = f0(index_t+1:Nft) -2*Nq_t;

Nz_save = Nz/2-1;
Nt_save = Nt/2-1;

sq_save = [1:Nt_save,Nft-Nt_save+2:Nft];

fts = ft(sq_save);

dft = ft(2)-ft(1);
dkz = kz(2)-kz(1);

%%% The length scale (It can be obtained from RANS, DNS, or experiments)
%%% The current file is from DNS. See Eq. (29) of Ref. (Ying, Li & Fu, JCP, 2025).
load('./Length_scale.mat','L')

%% E_iso (1D) from the spectra model

ftm = fts(1:Nt_save);
Nftm = length(ftm);

Umean1 = Umean;
Umean1(1) = 1e-6;  % Avoid Inf value at the wall when calculating ft./U (Specific value does not influence the result)

K =  (UU([2 2:Ny])+VV([2 2:Ny])+WW([2 2:Ny]))/2;
E_iso = f_Spectra_Pope(K',L',nu , fts(1:Nt_save)./Umean1');  % Nkx*Nkx*Nkx * Ny

%% Interpolate the E_HAT (3D) from E_iso (1D)
E_HAT = zeros(Nftm,Nkz,Nz_save,1,1,Ny);
for iiy = 2:Ny
    
    ratioU = ratio_X;
    kk3 = max(sqrt((ftm/Umean1(iiy)).^2 + (kz./ratioU)'.^2 + permute((kz(1:Nz_save)./ratioU),[3 2 1]).^2),1e-6);
    kk3_1d = kk3(:);
    
    E_HAT(:,:,:,1,1,iiy) = reshape(interp1(ftm/Umean1(iiy),E_iso(:,iiy),kk3_1d) , [Nftm,Nkz,Nz_save]);
    disp(['Assemble the 3-D spectra -Step-2... ',num2str(iiy),'/',num2str(Ny),'.'])
end
E_HAT(:,:,1,:,:,:) = 0.5*E_HAT(:,:,1,:,:,:);
E_HAT(find(isnan(E_HAT))) = 0;

E_HAT = E_HAT ./ sum(E_HAT*dft*dkz^2,[1 2 3]);

%% Scale E_HAT to the null space of wavenumbers
m_scale_C = zeros(Nftm,Nkz,Nz_save,3,3,Ny);
for iiy = 1:Ny
    
    scale_C = zeros(Nftm,Nkz,Nz_save,3,3);
    k1 = ftm/Umean1(iiy);
    k2 = kz';
    k3 = permute( kz(1:Nz_save) ,[3 2 1]);
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
    
    m_scale_C(:,:,:,:,:,iiy) = scale_C;

    disp(['Calculating the 3D spectra for UVW...',num2str(iiy),'/',num2str(Ny),'.']);

end

Cs_HAT = real(m_scale_C .* E_HAT(:,:,1:Nz_save,:,:,:)); % Nkx*Nkx*Nkx * 3 * 3 * Ny

%% Normalize the spectrum
uu = squeeze(1./mean(Cs_HAT(:,:,:,1,1,:),[1 2 3]));
vv = squeeze(1./mean(Cs_HAT(:,:,:,2,2,:),[1 2 3]));
ww = squeeze(1./mean(Cs_HAT(:,:,:,3,3,:),[1 2 3]));
uvw_sc = sqrt([uu,vv,ww]);
uvw_sc3 = permute(uvw_sc,[6 5 4 3 2 1]) .* permute(uvw_sc,[6 5 4 2 3 1]);
Cs_HAT_norm = Cs_HAT .* uvw_sc3;
Cs_HAT_norm(:,:,:,:,:,1) = Cs_HAT_norm(:,:,:,:,:,2);

%% Convert to physical space (y)

CST_H3 = permute(Cs_HAT_norm(:,:,1:Nz_save,:,:,:),[2 3 1 4 5 6]);% Y-Z-X -3-3 -Ny

CST_P1 = zeros(3*Ny,3*Ny,Nz_save,Nftm);

y3 = permute(y,[3 2 1]);
y4 = permute(y,[4 3 2 1]);
for iix = 1:Nftm
    
    m_fft = exp(1i*kz'.*(y3-y4))/Nkz;
    
    temp_3 = zeros(3*Ny,3*Ny,Nz_save);
    
    for ii1 = 1:3
        for ii2 = 1:3
            temp_1 = sqrt(squeeze(CST_H3(:,:,iix,ii1,ii2,:)) ...
                .* permute(squeeze(CST_H3(:,:,iix,ii1,ii2,:)),[1 2 4 3]));
            temp_2 = pagemtimes(m_fft,temp_1);
            
            temp_3( (ii1-1)*Ny+(1:Ny) , (ii2-1)*Ny+(1:Ny) , :) = permute(temp_2,[3 4 2 1]);
            
        end
    end
    
    CST_P1(:,:,:,iix) = temp_3;
    
    disp(['Constructing the CSD ... ',num2str(iix),'/',num2str(Nftm),'.'])
    
end

CST_P1 = (CST_P1 + conj(permute(CST_P1,[2 1 3 4])))/2;

CST_P2 = CST_P1;
CST_P2(1:Ny,Ny+(1:2*Ny),:,:) = 0;
CST_P2(Ny+(1:2*Ny),1:Ny,:,:) = 0;
CST_P2(Ny+(1:Ny),2*Ny+(1:Ny),:,:) = 0;
CST_P2(2*Ny+(1:Ny),Ny+(1:Ny),:,:) = 0;

%% Match the target Reynolds stress
Rij = zeros(3*Ny,3*Ny);
for iiy1 = 2:Ny
    Ruvw = [UU(iiy1) UV(iiy1)   0 ;
            UV(iiy1)   VV(iiy1) 0 ;
            0          0          WW(iiy1)];
    RL = chol(Ruvw);
    Rij(((1:3)-1)*Ny+iiy1,((1:3)-1)*Ny+iiy1) = RL;
end

CSD_temp1 = pagemtimes(Rij',CST_P2);
CSD_temp2 = pagemtimes(CSD_temp1,Rij);

% Rescale the Reynods stress profile
Euu_rec = zeros(Ny,1);
Evv_rec = zeros(Ny,1);
Euv_rec = zeros(Ny,1);
Eww_rec = zeros(Ny,1);
for ii1 = 1:Ny
    Euu_rec(ii1) = sum(CSD_temp2(ii1,ii1,:,:)*dft*dkz,[3 4]);
    Evv_rec(ii1) = sum(CSD_temp2(Ny+ii1,Ny+ii1,:,:)*dft*dkz,[3 4]);
    Euv_rec(ii1) = sum(CSD_temp2(Ny+ii1,ii1,:,:)*dft*dkz,[3 4]);
    Eww_rec(ii1) = sum(CSD_temp2(2*Ny+ii1,2*Ny+ii1,:,:)*dft*dkz,[3 4]);
end
factorR = max(UU)/max(Euu_rec);

Euu_rec = factorR*Euu_rec;
Evv_rec = factorR*Evv_rec;
Euv_rec = factorR*Euv_rec;
Eww_rec = factorR*Eww_rec;
CST = real(factorR*CSD_temp2);

%% Write the CSD tensor for GSI
save(['./CST_rx',num2str(ratio_X),'.mat'],'CST','-v7.3');

