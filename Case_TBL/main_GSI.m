% Reference paper: Ying, Li & Fu, Generalized synthetic inflow generation %
% method for divergence-free inhomogeneous turbulence, JCP, 2025          %
% Authors: Anjia Ying, Zhigang Li, Lin Fu                                 %

%%%%% Example program: Construct the eigenmodes with GSI, GSI-P and SFG for
%%%%% CFD computation of TBL
%%%%% Prerequisite: 3D covariance matrix stored in "CST_rx1.mat" obtained 
%%%%% by running "main_TRL.m" beforehand

clear;close all;

addpath('../Lib_GSI');

dir_ModeData = './Mode/';    % Direction where the divergence-free eigenmodes are stored
mkdir(dir_ModeData)


%%%% Set the number of workers used for parallel computation.
%%%% If not specified, the program will run with default settings.
% n_workers = 32;
% parpool(n_workers);

%% Mean quantities of velocities
load('./Mean_U_V.mat','Umean','Vmean','UU','VV','WW','UV')
Uinf = Umean(end);

U = Umean;
U(1) = 1e-6;  % Avoid Inf value at the wall when calculating ft./U (Specific value does not influence the result)

%% Grid settings (consistent with the CFD computation, as well as the main_TRL.m)
Nt = 400;
Ny = 151;
Nz = 150;

Lx = 200;
Lz = 15;

t = (1:Nt)/Nt*Lx / Uinf;
z = (1:Nz)/Nz*Lz;

y = load('./yp.dat');

d_y = (y(3:end)-y(1:end-2))/2;
d_y = [d_y(1);d_y;d_y(end)];

d_t = t(2)-t(1);
d_z = z(2)-z(1);

%% Wavenumbers and frequencies
Nkz = Nz;
Nft = Nt;

% Nyquist frequency
Nq_z = 2*pi/d_z/2;
Nq_t = 2*pi/d_t/2;

% Frequencies from 0 to 2Nq
kz0= 2*Nq_z*(0:Nkz-1)'/Nkz;
f0 = 2*Nq_t*(0:Nft-1)'/Nft;

% Restore the frequency ranging from -Nq to Nq
index_z = find(kz0 > Nq_z,1)-1;
index_t = find(f0  > Nq_t,1)-1;

kz = kz0;
kz(index_z+1:Nkz) = kz0(index_z+1:Nkz)-2*Nq_z;

ft = f0;
ft(index_t+1:Nft)  = f0(index_t+1:Nft) -2*Nq_t;

%
dft = ft(2)-ft(1);
dkz = kz(2)-kz(1);

%% Calculate differencing matrix in y direction
% Boundary conditions (BCs): lower bound: no-slip wall; upper bound: free
% Input: Coordinates; Accuracy order for differencing; BC at the starting point; BC at the end point
% Output: Differencing matrix
Dy = f_Diff(y,4,'no-slip wall','free');

%% Load the covariance matrix CST from DNS statistics or TRL (constructed in "main_TRL.m")

load('./CST_rx1.mat','CST');

[Nzs,Nts] = size(CST,[3 4]);

R_up = CST(:,:,1:Nzs,1:Nts);

R_um = R_up;
R_um(2*Ny+1:end,:,:,:) = -R_um(2*Ny+1:end,:,:,:);
R_um(:,2*Ny+1:end,:,:) = -R_um(:,2*Ny+1:end,:,:);

%% Initialize the phases of eigenmodes
rng(1);
sq_cal = 1:Ny;  % Sequence of the leading eigen modes restored for calculation

ny_cal = length(sq_cal);

% Random phases for SFG and the initialization of GSI-P
phase_random_s = random('uniform',0,2*pi,Nts,Nzs,ny_cal);
phase_random_c = random('uniform',0,2*pi,Nts,Nzs,ny_cal);

%% Initialize the output results
% Velocity
uu_GSI_P = zeros(Ny,Nz);
vv_GSI_P = zeros(Ny,Nz);
ww_GSI_P = zeros(Ny,Nz);

uu_GSI = zeros(Ny,Nz);
vv_GSI = zeros(Ny,Nz);
ww_GSI = zeros(Ny,Nz);

uu_SFG = zeros(Ny,Nz);
vv_SFG = zeros(Ny,Nz);
ww_SFG = zeros(Ny,Nz);

% Divergence
duu_GSI_P = zeros(Ny,Nz);
dvv_GSI_P = zeros(Ny,Nz);
dww_GSI_P = zeros(Ny,Nz);

duu_GSI = zeros(Ny,Nz);
dvv_GSI = zeros(Ny,Nz);
dww_GSI = zeros(Ny,Nz);

duu_SFG = zeros(Ny,Nz);
dvv_SFG = zeros(Ny,Nz);
dww_SFG = zeros(Ny,Nz);

% Energy spectra
E_u = zeros(Ny,1);
E_v = zeros(Ny,1);
E_w = zeros(Ny,1);

S_u_GSI_P = zeros(Ny,Nts,Nzs);
S_v_GSI_P = zeros(Ny,Nts,Nzs);
S_w_GSI_P = zeros(Ny,Nts,Nzs);

S_u_GSI = zeros(Ny,Nts,Nzs);
S_v_GSI = zeros(Ny,Nts,Nzs);
S_w_GSI = zeros(Ny,Nts,Nzs);

S_u_SFG = zeros(Ny,Nts,Nzs);
S_v_SFG = zeros(Ny,Nts,Nzs);
S_w_SFG = zeros(Ny,Nts,Nzs);

it_test = 1; % The time index to display the sample instantaneous field, which does not influence the output eigenmodes.

%% Construct the divergence-free eigenmodes
for iix = 1:Nts
    
    % Temporal variables for calculating instantaneous velocity
    sx = sin(ft(iix)*t(it_test));
    cx = cos(ft(iix)*t(it_test));
    
    % Temporal variables for sampling the spectral energy
    E_ut = zeros(Ny,1);
    E_vt = zeros(Ny,1);
    E_wt = zeros(Ny,1);
    
    S_u_GSIt_P = zeros(Ny,1,Nzs);
    S_v_GSIt_P = zeros(Ny,1,Nzs);
    S_w_GSIt_P = zeros(Ny,1,Nzs);

    S_u_GSIt = zeros(Ny,1,Nzs);
    S_v_GSIt = zeros(Ny,1,Nzs);
    S_w_GSIt = zeros(Ny,1,Nzs);

    S_u_SFGt = zeros(Ny,1,Nzs);
    S_v_SFGt = zeros(Ny,1,Nzs);
    S_w_SFGt = zeros(Ny,1,Nzs);
    
    %% Initialize the eigenmodes for each frequency
    eig_Us_GSI_P = zeros(Ny,Nz);
    eig_Uc_GSI_P = zeros(Ny,Nz);
    eig_Vs_GSI_P = zeros(Ny,Nz);
    eig_Vc_GSI_P = zeros(Ny,Nz);
    eig_Ws_GSI_P = zeros(Ny,Nz);
    eig_Wc_GSI_P = zeros(Ny,Nz);
    
    %
    eig_Us_GSI = zeros(Ny,Nz);
    eig_Uc_GSI = zeros(Ny,Nz);
    eig_Vs_GSI = zeros(Ny,Nz);
    eig_Vc_GSI = zeros(Ny,Nz);
    eig_Ws_GSI = zeros(Ny,Nz);
    eig_Wc_GSI = zeros(Ny,Nz);

    eig_Us_SFG = zeros(Ny,Nz);
    eig_Uc_SFG = zeros(Ny,Nz);
    eig_Vs_SFG = zeros(Ny,Nz);
    eig_Vc_SFG = zeros(Ny,Nz);
    eig_Ws_SFG = zeros(Ny,Nz);
    eig_Wc_SFG = zeros(Ny,Nz);
    
    parfor iiz = 1:Nzs
        
        % Skip the zero-wavenumber case corresponding to the mean state
        if (ft(iix) == 0 && kz(iiz) == 0)
            continue
        end

        sz = sin(kz(iiz)*z);
        cz = cos(kz(iiz)*z);
        s_s = sx.*sz;
        s_c = sx.*cz;
        c_s = cx.*sz;
        c_c = cx.*cz;

        %%%%% Regularization term for eliminating the divergence %%%%%
        gamma_U = 200;     % Regularization coefficient
        
        Ez = kz(iiz)*eye(Ny);
        Ex = diag( -ft(iix)./U )*eye(Ny);
        % ss
        ExR = Ex;
        EzR = Ez;
        D3_SS = [ExR , Dy , EzR].*sqrt(d_y);
        IE3 = (eye(3*Ny) + gamma_U*D3_SS'*D3_SS);
        inv_GSI_SS = IE3\eye(3*Ny);
        % sc
        ExR = Ex;
        EzR = -Ez;
        D3_SC = [ExR , Dy , EzR].*sqrt(d_y);
        IE3 = (eye(3*Ny) + gamma_U*D3_SC'*D3_SC);
        inv_GSI_SC = IE3\eye(3*Ny);
        % cs
        ExR = -Ex;
        EzR = Ez;
        D3_CS = [ExR , Dy , EzR].*sqrt(d_y);
        IE3 = (eye(3*Ny) + gamma_U*D3_CS'*D3_CS);
        inv_GSI_CS = IE3\eye(3*Ny);
        % cc
        ExR = -Ex;
        EzR = -Ez;
        D3_CC = [ExR , Dy , EzR].*sqrt(d_y);
        IE3 = (eye(3*Ny) + gamma_U*D3_CC'*D3_CC);
        inv_GSI_CC = IE3\eye(3*Ny);

        %%% The bases for synthesizing the inflow turbulence
        % Part 1
        CSD_u3 = R_up(:,:,iiz,iix)*dft*dkz;

        dy3 = [d_y;d_y;d_y];

        [eig_U0,eig_SU0] = eig(CSD_u3.*sqrt(dy3*dy3'));

        eig_SU = diag(abs(eig_SU0));    eig_U1 = (eig_U0./sqrt(dy3)).*sqrt(eig_SU)';

        [~,sq_eigU] = sort(eig_SU); eig_U2 = eig_U1(:,flipud(sq_eigU));

        eig_Up = eig_U2(1:Ny,sq_cal);
        eig_Vp = eig_U2(Ny+(1:Ny),sq_cal);
        eig_Wp = eig_U2(2*Ny+(1:Ny),sq_cal);
        
        % Part2
        CSD_u3 = R_um(:,:,iiz,iix)*dft*dkz;

        dy3 = [d_y;d_y;d_y];

        [eig_U0,eig_SU0] = eig(CSD_u3.*sqrt(dy3*dy3'));

        eig_SU = diag(abs(eig_SU0));    eig_U1 = (eig_U0./sqrt(dy3)).*sqrt(eig_SU)';

        [~,sq_eigU] = sort(eig_SU); eig_U2 = eig_U1(:,flipud(sq_eigU));

        eig_Um = eig_U2(1:Ny,sq_cal);
        eig_Vm = eig_U2(Ny+(1:Ny),sq_cal);
        eig_Wm = eig_U2(2*Ny+(1:Ny),sq_cal);
        
        %%% Target values (for validating the generation results)
        Eut = real(diag(CSD_u3(1:Ny,1:Ny)));
        Evt = real(diag(CSD_u3(Ny+(1:Ny),Ny+(1:Ny))));
        Ewt = real(diag(CSD_u3(2*Ny+(1:Ny),2*Ny+(1:Ny))));
        
        E_ut = E_ut + Eut;
        E_vt = E_vt + Evt;
        E_wt = E_wt + Ewt;
        
        %% GSI
        aRp = real(eig_Up);
        rRp = real(eig_Wp);
        bRp = real(eig_Vp);
        aIp = imag(eig_Up);
        rIp = imag(eig_Wp);
        bIp = imag(eig_Vp);

        aRm = real(eig_Um);
        rRm = real(eig_Wm);
        bRm = real(eig_Vm);
        aIm = imag(eig_Um);
        rIm = imag(eig_Wm);
        bIm = imag(eig_Vm);
        

        %%%%% GSI step 1: Phase optimization %%%%%

        %%% Optimize the phases
        phase_init = [ squeeze(phase_random_c(iix,iiz,sq_cal)) ; squeeze(phase_random_s(iix,iiz,sq_cal)) ]' ;
        [phase_opt,J,J_init,~] = f_GSI_1D(eig_Up,eig_Vp,eig_Wp,eig_Um,eig_Vm,eig_Wm,ft(iix),kz(iiz),U,Dy,d_y,phase_init(:));

        %%% Assemble the eigenmodes with optimized phases
        phase_xp = phase_opt(1,:);
        phase_yp = phase_opt(1,:);
        phase_zp = phase_opt(1,:);

        phase_xm = phase_opt(2,:);
        phase_ym = phase_opt(2,:);
        phase_zm = phase_opt(2,:);

        c_xp = cos(phase_xp);
        s_xp = sin(phase_xp);
        c_yp = cos(phase_yp);
        s_yp = sin(phase_yp);
        c_zp = cos(phase_zp);
        s_zp = sin(phase_zp);

        c_xm = cos(phase_xm);
        s_xm = sin(phase_xm);
        c_ym = cos(phase_ym);
        s_ym = sin(phase_ym);
        c_zm = cos(phase_zm);
        s_zm = sin(phase_zm);

        uu_cc = permute(sum( ( aRp.*( c_xp)-aIp.*( s_xp) + aRm.*( c_xm)-aIm.*( s_xm) ) , 2),[1 3 2]);
        uu_ss = permute(sum( ( aRp.*(-c_xp)-aIp.*(-s_xp) + aRm.*( c_xm)-aIm.*( s_xm) ) , 2),[1 3 2]);
        uu_sc = permute(sum( ( aRp.*(-s_xp)-aIp.*( c_xp) + aRm.*(-s_xm)-aIm.*( c_xm) ) , 2),[1 3 2]);
        uu_cs = permute(sum( ( aRp.*(-s_xp)-aIp.*( c_xp) + aRm.*( s_xm)-aIm.*(-c_xm) ) , 2),[1 3 2]);

        vv_cc = permute(sum( ( bRp.*( c_yp)-bIp.*( s_yp) + bRm.*( c_ym)-bIm.*( s_ym) ) , 2),[1 3 2]);
        vv_ss = permute(sum( ( bRp.*(-c_yp)-bIp.*(-s_yp) + bRm.*( c_ym)-bIm.*( s_ym) ) , 2),[1 3 2]);
        vv_sc = permute(sum( ( bRp.*(-s_yp)-bIp.*( c_yp) + bRm.*(-s_ym)-bIm.*( c_ym) ) , 2),[1 3 2]);
        vv_cs = permute(sum( ( bRp.*(-s_yp)-bIp.*( c_yp) + bRm.*( s_ym)-bIm.*(-c_ym) ) , 2),[1 3 2]);

        ww_cc = permute(sum( ( rRp.*( c_zp)-rIp.*( s_zp) + rRm.*( c_zm)-rIm.*( s_zm) ) , 2),[1 3 2]);
        ww_ss = permute(sum( ( rRp.*(-c_zp)-rIp.*(-s_zp) + rRm.*( c_zm)-rIm.*( s_zm) ) , 2),[1 3 2]);
        ww_sc = permute(sum( ( rRp.*(-s_zp)-rIp.*( c_zp) + rRm.*(-s_zm)-rIm.*( c_zm) ) , 2),[1 3 2]);
        ww_cs = permute(sum( ( rRp.*(-s_zp)-rIp.*( c_zp) + rRm.*( s_zm)-rIm.*(-c_zm) ) , 2),[1 3 2]);

        % 
        eig_Us_GSI_P = eig_Us_GSI_P + (sz.*uu_ss + cz.*uu_sc);
        eig_Uc_GSI_P = eig_Uc_GSI_P + (sz.*uu_cs + cz.*uu_cc);
        eig_Vs_GSI_P = eig_Vs_GSI_P + (sz.*vv_ss + cz.*vv_sc);
        eig_Vc_GSI_P = eig_Vc_GSI_P + (sz.*vv_cs + cz.*vv_cc);
        eig_Ws_GSI_P = eig_Ws_GSI_P + (sz.*ww_ss + cz.*ww_sc);
        eig_Wc_GSI_P = eig_Wc_GSI_P + (sz.*ww_cs + cz.*ww_cc);
        
        % Divergence (for validating the generation results)
        duu_GSI_P = duu_GSI_P + (-ft(iix)./U) .* (c_s.*uu_ss + c_c.*uu_sc - s_s.*uu_cs - s_c.*uu_cc);
        dvv_GSI_P = dvv_GSI_P + Dy*(s_s.*vv_ss + s_c.*vv_sc + c_s.*vv_cs + c_c.*vv_cc);
        dww_GSI_P = dww_GSI_P + ( kz(iiz)   ) .* (s_c.*ww_ss - s_s.*ww_sc + c_c.*ww_cs - c_s.*ww_cc);
        
        % Energy spectra (for validating the generation results)
        S_u_GSIt_P(:,1,iiz) = (uu_ss.^2 + uu_sc.^2 + uu_cs.^2 + uu_cc.^2)/4;
        S_v_GSIt_P(:,1,iiz) = (vv_ss.^2 + vv_sc.^2 + vv_cs.^2 + vv_cc.^2)/4;
        S_w_GSIt_P(:,1,iiz) = (ww_ss.^2 + ww_sc.^2 + ww_cs.^2 + ww_cc.^2)/4;

        %%%%% GSI step2: Velocity adjustment %%%%%

        % Adjust the eigenmodes to eliminate divergence
        uvw_ss = inv_GSI_CC*([uu_cs(:);vv_ss(:);ww_sc(:)]);
        uvw_sc = inv_GSI_CS*([uu_cc(:);vv_sc(:);ww_ss(:)]);
        uvw_cs = inv_GSI_SC*([uu_ss(:);vv_cs(:);ww_cc(:)]);
        uvw_cc = inv_GSI_SS*([uu_sc(:);vv_cc(:);ww_cs(:)]);

        uu2_cs = uvw_ss(1:Ny);vv2_ss = uvw_ss(Ny+1:2*Ny);ww2_sc = uvw_ss(2*Ny+1:end);
        uu2_cc = uvw_sc(1:Ny);vv2_sc = uvw_sc(Ny+1:2*Ny);ww2_ss = uvw_sc(2*Ny+1:end);
        uu2_ss = uvw_cs(1:Ny);vv2_cs = uvw_cs(Ny+1:2*Ny);ww2_cc = uvw_cs(2*Ny+1:end);
        uu2_sc = uvw_cc(1:Ny);vv2_cc = uvw_cc(Ny+1:2*Ny);ww2_cs = uvw_cc(2*Ny+1:end);

        eig_Us_GSI = eig_Us_GSI + (sz.*uu2_ss + cz.*uu2_sc);
        eig_Uc_GSI = eig_Uc_GSI + (sz.*uu2_cs + cz.*uu2_cc);
        eig_Vs_GSI = eig_Vs_GSI + (sz.*vv2_ss + cz.*vv2_sc);
        eig_Vc_GSI = eig_Vc_GSI + (sz.*vv2_cs + cz.*vv2_cc);
        eig_Ws_GSI = eig_Ws_GSI + (sz.*ww2_ss + cz.*ww2_sc);
        eig_Wc_GSI = eig_Wc_GSI + (sz.*ww2_cs + cz.*ww2_cc);
        
        % Divergence (for validating the generation results)
        duu_GSI = duu_GSI + (-ft(iix)./U) .* (c_s.*uu2_ss + c_c.*uu2_sc - s_s.*uu2_cs - s_c.*uu2_cc);
        dvv_GSI = dvv_GSI + Dy*(s_s.*vv2_ss + s_c.*vv2_sc + c_s.*vv2_cs + c_c.*vv2_cc);
        dww_GSI = dww_GSI + ( kz(iiz)   ) .* (s_c.*ww2_ss - s_s.*ww2_sc + c_c.*ww2_cs - c_s.*ww2_cc);
        
        % Energy spectra (for validating the generation results)
        S_u_GSIt(:,1,iiz) = (uu2_ss.^2 + uu2_sc.^2 + uu2_cs.^2 + uu2_cc.^2)/4;
        S_v_GSIt(:,1,iiz) = (vv2_ss.^2 + vv2_sc.^2 + vv2_cs.^2 + vv2_cc.^2)/4;
        S_w_GSIt(:,1,iiz) = (ww2_ss.^2 + ww2_sc.^2 + ww2_cs.^2 + ww2_cc.^2)/4;
        

        %% Reference case: SFG (without divergence-eliminating operations)

        %%% Assemble the eigenmodes with random phases
        phase_xp = permute(phase_random_c(iix,iiz,:),[1 3 2]);
        phase_yp = permute(phase_random_c(iix,iiz,:),[1 3 2]);
        phase_zp = permute(phase_random_c(iix,iiz,:),[1 3 2]);

        phase_xm = permute(phase_random_s(iix,iiz,:),[1 3 2]);
        phase_ym = permute(phase_random_s(iix,iiz,:),[1 3 2]);
        phase_zm = permute(phase_random_s(iix,iiz,:),[1 3 2]);

        %
        c_xp = cos(phase_xp);
        s_xp = sin(phase_xp);
        c_yp = cos(phase_yp);
        s_yp = sin(phase_yp);
        c_zp = cos(phase_zp);
        s_zp = sin(phase_zp);

        c_xm = cos(phase_xm);
        s_xm = sin(phase_xm);
        c_ym = cos(phase_ym);
        s_ym = sin(phase_ym);
        c_zm = cos(phase_zm);
        s_zm = sin(phase_zm);

        uuR_cc = permute(sum( ( aRp.*( c_xp)-aIp.*( s_xp) + aRm.*( c_xm)-aIm.*( s_xm) ) , 2),[1 3 2]);
        uuR_ss = permute(sum( ( aRp.*(-c_xp)-aIp.*(-s_xp) + aRm.*( c_xm)-aIm.*( s_xm) ) , 2),[1 3 2]);
        uuR_sc = permute(sum( ( aRp.*(-s_xp)-aIp.*( c_xp) + aRm.*(-s_xm)-aIm.*( c_xm) ) , 2),[1 3 2]);
        uuR_cs = permute(sum( ( aRp.*(-s_xp)-aIp.*( c_xp) + aRm.*( s_xm)-aIm.*(-c_xm) ) , 2),[1 3 2]);

        vvR_cc = permute(sum( ( bRp.*( c_yp)-bIp.*( s_yp) + bRm.*( c_ym)-bIm.*( s_ym) ) , 2),[1 3 2]);
        vvR_ss = permute(sum( ( bRp.*(-c_yp)-bIp.*(-s_yp) + bRm.*( c_ym)-bIm.*( s_ym) ) , 2),[1 3 2]);
        vvR_sc = permute(sum( ( bRp.*(-s_yp)-bIp.*( c_yp) + bRm.*(-s_ym)-bIm.*( c_ym) ) , 2),[1 3 2]);
        vvR_cs = permute(sum( ( bRp.*(-s_yp)-bIp.*( c_yp) + bRm.*( s_ym)-bIm.*(-c_ym) ) , 2),[1 3 2]);

        wwR_cc = permute(sum( ( rRp.*( c_zp)-rIp.*( s_zp) + rRm.*( c_zm)-rIm.*( s_zm) ) , 2),[1 3 2]);
        wwR_ss = permute(sum( ( rRp.*(-c_zp)-rIp.*(-s_zp) + rRm.*( c_zm)-rIm.*( s_zm) ) , 2),[1 3 2]);
        wwR_sc = permute(sum( ( rRp.*(-s_zp)-rIp.*( c_zp) + rRm.*(-s_zm)-rIm.*( c_zm) ) , 2),[1 3 2]);
        wwR_cs = permute(sum( ( rRp.*(-s_zp)-rIp.*( c_zp) + rRm.*( s_zm)-rIm.*(-c_zm) ) , 2),[1 3 2]);

        %
        eig_Us_SFG = eig_Us_SFG + (sz.*uuR_ss + cz.*uuR_sc);
        eig_Uc_SFG = eig_Uc_SFG + (sz.*uuR_cs + cz.*uuR_cc);
        eig_Vs_SFG = eig_Vs_SFG + (sz.*vvR_ss + cz.*vvR_sc);
        eig_Vc_SFG = eig_Vc_SFG + (sz.*vvR_cs + cz.*vvR_cc);
        eig_Ws_SFG = eig_Ws_SFG + (sz.*wwR_ss + cz.*wwR_sc);
        eig_Wc_SFG = eig_Wc_SFG + (sz.*wwR_cs + cz.*wwR_cc);
        
        % Divergence (for validating the generation results)
        duu_SFG = duu_SFG + (-ft(iix)./U) .* (c_s.*uuR_ss + c_c.*uuR_sc - s_s.*uuR_cs - s_c.*uuR_cc);
        dvv_SFG = dvv_SFG + Dy*(s_s.*vvR_ss + s_c.*vvR_sc + c_s.*vvR_cs + c_c.*vvR_cc);
        dww_SFG = dww_SFG + ( kz(iiz)   ) .* (s_c.*wwR_ss - s_s.*wwR_sc + c_c.*wwR_cs - c_s.*wwR_cc);
        
        % Energy spectra (for validating the generation results)
        S_u_SFGt(:,1,iiz) = (uuR_ss.^2 + uuR_sc.^2 + uuR_cs.^2 + uuR_cc.^2)/4;
        S_v_SFGt(:,1,iiz) = (vvR_ss.^2 + vvR_sc.^2 + vvR_cs.^2 + vvR_cc.^2)/4;
        S_w_SFGt(:,1,iiz) = (wwR_ss.^2 + wwR_sc.^2 + wwR_cs.^2 + wwR_cc.^2)/4;

        disp(['Calculating the parameters -> X: ',num2str(iix),', Z: ',num2str(iiz),'.']);
        
    end

    uu_GSI_P = uu_GSI_P + (sx.*eig_Us_GSI_P + cx.*eig_Uc_GSI_P);
    vv_GSI_P = vv_GSI_P + (sx.*eig_Vs_GSI_P + cx.*eig_Vc_GSI_P);
    ww_GSI_P = ww_GSI_P + (sx.*eig_Ws_GSI_P + cx.*eig_Wc_GSI_P);

    uu_GSI = uu_GSI + (sx.*eig_Us_GSI + cx.*eig_Uc_GSI);
    vv_GSI = vv_GSI + (sx.*eig_Vs_GSI + cx.*eig_Vc_GSI);
    ww_GSI = ww_GSI + (sx.*eig_Ws_GSI + cx.*eig_Wc_GSI);

    uu_SFG = uu_SFG + (sx.*eig_Us_SFG + cx.*eig_Uc_SFG);
    vv_SFG = vv_SFG + (sx.*eig_Vs_SFG + cx.*eig_Vc_SFG);
    ww_SFG = ww_SFG + (sx.*eig_Ws_SFG + cx.*eig_Wc_SFG);

    E_u = E_u + E_ut;
    E_v = E_v + E_vt;
    E_w = E_w + E_wt;

    S_u_GSI_P(:,iix,:) = S_u_GSIt_P;
    S_v_GSI_P(:,iix,:) = S_v_GSIt_P;
    S_w_GSI_P(:,iix,:) = S_w_GSIt_P;

    S_u_GSI(:,iix,:) = S_u_GSIt;
    S_v_GSI(:,iix,:) = S_v_GSIt;
    S_w_GSI(:,iix,:) = S_w_GSIt;

    S_u_SFG(:,iix,:) = S_u_SFGt;
    S_v_SFG(:,iix,:) = S_v_SFGt;
    S_w_SFG(:,iix,:) = S_w_SFGt;
    
    %% Save the GSI data
    % GSI-P
    fid=fopen([dir_ModeData,'uuS_GSI_P_iX-',num2str(iix),'.bin'],'w');
    fwrite(fid,eig_Us_GSI_P(:),'float32');
    fclose(fid);

    fid=fopen([dir_ModeData,'uuC_GSI_P_iX-',num2str(iix),'.bin'],'w');
    fwrite(fid,eig_Uc_GSI_P(:),'float32');
    fclose(fid);

    fid=fopen([dir_ModeData,'vvS_GSI_P_iX-',num2str(iix),'.bin'],'w');
    fwrite(fid,eig_Vs_GSI_P(:),'float32');
    fclose(fid);

    fid=fopen([dir_ModeData,'vvC_GSI_P_iX-',num2str(iix),'.bin'],'w');
    fwrite(fid,eig_Vc_GSI_P(:),'float32');
    fclose(fid);

    fid=fopen([dir_ModeData,'wwS_GSI_P_iX-',num2str(iix),'.bin'],'w');
    fwrite(fid,eig_Ws_GSI_P(:),'float32');
    fclose(fid);

    fid=fopen([dir_ModeData,'wwC_GSI_P_iX-',num2str(iix),'.bin'],'w');
    fwrite(fid,eig_Wc_GSI_P(:),'float32');
    fclose(fid);
    
    % GSI-P&V (GSI for brevity)
    fid=fopen([dir_ModeData,'uuS_GSI_iX-',num2str(iix),'.bin'],'w');
    fwrite(fid,eig_Us_GSI(:),'float32');
    fclose(fid);

    fid=fopen([dir_ModeData,'uuC_GSI_iX-',num2str(iix),'.bin'],'w');
    fwrite(fid,eig_Uc_GSI(:),'float32');
    fclose(fid);

    fid=fopen([dir_ModeData,'vvS_GSI_iX-',num2str(iix),'.bin'],'w');
    fwrite(fid,eig_Vs_GSI(:),'float32');
    fclose(fid);

    fid=fopen([dir_ModeData,'vvC_GSI_iX-',num2str(iix),'.bin'],'w');
    fwrite(fid,eig_Vc_GSI(:),'float32');
    fclose(fid);

    fid=fopen([dir_ModeData,'wwS_GSI_iX-',num2str(iix),'.bin'],'w');
    fwrite(fid,eig_Ws_GSI(:),'float32');
    fclose(fid);

    fid=fopen([dir_ModeData,'wwC_GSI_iX-',num2str(iix),'.bin'],'w');
    fwrite(fid,eig_Wc_GSI(:),'float32');
    fclose(fid);
    
    % SFG
    fid=fopen([dir_ModeData,'uuS_SFG_iX-',num2str(iix),'.bin'],'w');
    fwrite(fid,eig_Us_SFG(:),'float32');
    fclose(fid);

    fid=fopen([dir_ModeData,'uuC_SFG_iX-',num2str(iix),'.bin'],'w');
    fwrite(fid,eig_Uc_SFG(:),'float32');
    fclose(fid);

    fid=fopen([dir_ModeData,'vvS_SFG_iX-',num2str(iix),'.bin'],'w');
    fwrite(fid,eig_Vs_SFG(:),'float32');
    fclose(fid);

    fid=fopen([dir_ModeData,'vvC_SFG_iX-',num2str(iix),'.bin'],'w');
    fwrite(fid,eig_Vc_SFG(:),'float32');
    fclose(fid);

    fid=fopen([dir_ModeData,'wwS_SFG_iX-',num2str(iix),'.bin'],'w');
    fwrite(fid,eig_Ws_SFG(:),'float32');
    fclose(fid);

    fid=fopen([dir_ModeData,'wwC_SFG_iX-',num2str(iix),'.bin'],'w');
    fwrite(fid,eig_Wc_SFG(:),'float32');
    fclose(fid);
    
end


%% Test the generated results

%%% The integrated variance of velocity divergence
div_GSI = sum(sum((duu_GSI+dvv_GSI+dww_GSI).^2.*d_y,1));
div_GSI_P = sum(sum((duu_GSI_P+dvv_GSI_P+dww_GSI_P).^2.*d_y,1));
div_SFG = sum(sum((duu_SFG+dvv_SFG+dww_SFG).^2.*d_y,1));

disp(['The divergence results are: GSI: ',num2str(div_GSI),'; GSI-P: ',num2str(div_GSI_P),'; SFG: ',num2str(div_SFG),'. DIV_GSI/DIV_SFG = ',num2str(div_GSI/div_SFG),'.'])

%%% The energy profile
E_u_GSI_P = sum(S_u_GSI_P,[2 3]);
E_v_GSI_P = sum(S_v_GSI_P,[2 3]);
E_w_GSI_P = sum(S_w_GSI_P,[2 3]);

E_u_GSI = sum(S_u_GSI,[2 3]);
E_v_GSI = sum(S_v_GSI,[2 3]);
E_w_GSI = sum(S_w_GSI,[2 3]);

E_u_SFG = sum(S_u_SFG,[2 3]);
E_v_SFG = sum(S_v_SFG,[2 3]);
E_w_SFG = sum(S_w_SFG,[2 3]);

figure('Position',[100, 100, 800, 600]);
E_all = {E_u, E_v, E_w};
E_GSI_P = {E_u_GSI_P, E_v_GSI_P, E_w_GSI_P};
E_GSI = {E_u_GSI, E_v_GSI, E_w_GSI};
E_SFG = {E_u_SFG, E_v_SFG, E_w_SFG};
xlabels = {'$\left\langle u^{\prime} u^{\prime} \right\rangle$', ...
           '$\left\langle v^{\prime} v^{\prime} \right\rangle$', ...
           '$\left\langle w^{\prime} w^{\prime} \right\rangle$'};
ylabels = {'$y$', '$y$', '$y$'};

for i = 1:3
    subplot(1,3,i)
    plot(sqrt(E_all{i}), y, 'k-'); hold on
    plot(sqrt(E_GSI{i}), y, 'r-');
    plot(sqrt(E_GSI_P{i}), y, 'g--');
    plot(sqrt(E_SFG{i}), y, 'b.');
    xlim([0 1.5*max(sqrt(E_GSI_P{i}))])
    xlabel(xlabels{i}, 'Interpreter', 'latex')
    ylabel(ylabels{i}, 'Interpreter', 'latex')
    if i == 1
        legend('Target','GSI-P\&V','GSI-P','SFG', 'Interpreter', 'latex')
    end
end

savefig('Energy_profile_TBL.fig')

%%% Instantaneous velocity field
[zz,yy] = meshgrid(z,y);

figure('Position', [200, 200, 1500, 600]);

data = {uu_GSI(:,:,1),   vv_GSI(:,:,1),   ww_GSI(:,:,1),   duu_GSI+dvv_GSI+dww_GSI, ...
        uu_GSI_P(:,:,1), vv_GSI_P(:,:,1), ww_GSI_P(:,:,1), duu_GSI_P+dvv_GSI_P+dww_GSI_P, ...
        uu_SFG(:,:,1),   vv_SFG(:,:,1),   ww_SFG(:,:,1),   duu_SFG+dvv_SFG+dww_SFG};

titles = {'$u$ (GSI-P\&V)', '$v$ (GSI-P\&V)', '$w$ (GSI-P\&V)', 'DIV (GSI-P\&V)', ...
          '$u$ (GSI-P)', '$v$ (GSI-P)', '$w$ (GSI-P)', 'DIV (GSI-P)', ...
          '$u$ (SFG)', '$v$ (SFG)', '$w$ (SFG)', 'DIV (SFG)'};

caxis_limits = {
    [-1 1]*max(abs(uu_GSI_P),[],'all')*0.8, ...
    [-1 1]*max(abs(vv_GSI_P),[],'all')*0.8, ...
    [-1 1]*max(abs(ww_GSI_P),[],'all')*0.8, ...
    [-1 1]*mean(abs(duu_SFG+dvv_SFG+dww_SFG),'all'), ...
    [-1 1]*max(abs(uu_GSI_P),[],'all')*0.8, ...
    [-1 1]*max(abs(vv_GSI_P),[],'all')*0.8, ...
    [-1 1]*max(abs(ww_GSI_P),[],'all')*0.8, ...
    [-1 1]*mean(abs(duu_SFG+dvv_SFG+dww_SFG),'all'), ...
    [-1 1]*max(abs(uu_SFG),[],'all')*0.8, ...
    [-1 1]*max(abs(vv_SFG),[],'all')*0.8, ...
    [-1 1]*max(abs(ww_SFG),[],'all')*0.8, ...
    [-1 1]*mean(abs(duu_SFG+dvv_SFG+dww_SFG),'all')
    };

for i = 1:12
    subplot(3,4,i)
    pcolor(zz,yy,data{i});
    shading interp;
    daspect([1 1 1]);
    ylim([0 10]);
    colorbar;
    caxis(caxis_limits{i});
    title(titles{i}, 'Interpreter', 'latex');
    xlabel('$z$', 'Interpreter', 'latex');
    ylabel('$y$', 'Interpreter', 'latex');
    colormap(jet)
end

savefig('Instantaneous_field_TBL.fig')
