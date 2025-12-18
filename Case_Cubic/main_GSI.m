% Reference paper: Ying, Li & Fu, Generalized synthetic inflow generation %
% method for divergence-free inhomogeneous turbulence, JCP, 2025          %
% Authors: Anjia Ying, Zhigang Li, Lin Fu                                 %

%%%%% Example program: Construct the eigenmodes with GSI, GSI-P and SFG for
%%%%% CFD computation of homogeneous turbulence (inhomogeneous boundaries can be present)

clear;close all;
addpath('../Lib_GSI');

dir_ModeData = './Mode/';    % Direction where the generated eigenmodes are stored
mkdir(dir_ModeData);

%%%% Set the number of workers used for parallel computation.
%%%% If not specified, the program will run with default settings.
% n_workers = 32;
% parpool(n_workers);

%% Flow quantities

U = 10;      % Mean streamwise velocity
nu = 1/1000; % Kinetic viscosity

iic = 2;    % Here the second case (Field 1 in Section 3) is chosen for instance

% Parameters to construct the target turbulence properties
sigma_U = [1,  1  ,  1  ,  1  ,  2/3,  2/3,  2/3];
sigma_R = [0,  0  , -1/2, -1/2,  0  , -1/2, -1/2];
sigma_L = [1,  2/3,  1  ,  2/3,  1  ,  1  ,  2/3];

UU = 1;
VV = 1*sigma_U(iic);
WW = 1*sigma_U(iic)^2;
UV   = sqrt(sigma_U(iic))*sigma_R(iic);

% Reynolds stress
Ruvw = [UU       UV       0 ;
        UV       VV       0 ;
        0        0        WW ];

% Length scales
L1 = 1;
L2 = 1*sigma_L(iic);
L3 = 1*sigma_L(iic)^2;

mL = [L1,L2,L3];

%% Grid settings (consistent with the CFD computation)
Nt = 300;
Nx = 300;
Ny = 200;
Nz = 200;

ny = Ny;
nz = Nz;

x = (1:Nx)'/Nx*(3*pi);
z = (1:Nz) /Nz*(2*pi);
y = (1:Ny)'/Ny*(2*pi);

d_x = x(2)-x(1);
d_y = y(2)-y(1);
d_t = d_x/U;

t = d_t*(1:Nt);

%% Reconstruct the 3D cross-spectra with TRL model
data = f_TRL(z,t,U,nu,Ruvw,mL);

%% Wavenumbers and frequencies
ft      = data.ft;
kz      = data.kz;
Nt_mod = data.Nt_mod;
Nz_mod = data.Nz_mod;
CST = data.CST/100;

EST = squeeze(CST(1,1,:,:,:));

%% Differencing matrix in y direction
% Name of boundaries: 'periodic', 'symmetry', 'no-slip wall', 'free'
Dy = f_Diff(y,4,'periodic','periodic');

%% Initialize the phases of eigenmodes
rng(1);
phase_random_s = random('uniform',0,2*pi,Nt_mod,Nz_mod,Nz_mod,6);
phase_random_c = random('uniform',0,2*pi,Nt_mod,Nz_mod,Nz_mod,6);

%% Initialize the output results
% Velocity signal
uu_GSI_P = zeros(ny,nz);
vv_GSI_P = zeros(ny,nz);
ww_GSI_P = zeros(ny,nz);

uu_GSI = zeros(ny,nz);
vv_GSI = zeros(ny,nz);
ww_GSI = zeros(ny,nz);

uu_SFG = zeros(ny,nz);
vv_SFG = zeros(ny,nz);
ww_SFG = zeros(ny,nz);

% Divergence
duu_GSI_P = zeros(ny,nz);
dvv_GSI_P = zeros(ny,nz);
dww_GSI_P = zeros(ny,nz);

duu_GSI = zeros(ny,nz);
dvv_GSI = zeros(ny,nz);
dww_GSI = zeros(ny,nz);

duu_SFG = zeros(ny,nz);
dvv_SFG = zeros(ny,nz);
dww_SFG = zeros(ny,nz);

% Energy spectra
S_u_GSI_P = zeros(Nt_mod,Nz_mod,Nz_mod);
S_v_GSI_P = zeros(Nt_mod,Nz_mod,Nz_mod);
S_w_GSI_P = zeros(Nt_mod,Nz_mod,Nz_mod);

S_u_GSI = zeros(Nt_mod,Nz_mod,Nz_mod);
S_v_GSI = zeros(Nt_mod,Nz_mod,Nz_mod);
S_w_GSI = zeros(Nt_mod,Nz_mod,Nz_mod);

S_u_SFG = zeros(Nt_mod,Nz_mod,Nz_mod);
S_v_SFG = zeros(Nt_mod,Nz_mod,Nz_mod);
S_w_SFG = zeros(Nt_mod,Nz_mod,Nz_mod);

it_test = 1; % The time index to display the sample instantaneous field, which does not influence the output eigenmodes.

%% Construct the divergence-free eigenmodes
for iix = 1:Nt_mod
    
    % Temporal variables for calculating instantaneous velocity
    sx = sin(ft(iix)*t(it_test));
    cx = cos(ft(iix)*t(it_test));
    
    %% Initialize the eigenmodes for each frequency
    eig_Us_GSI_P = zeros(ny,nz);
    eig_Uc_GSI_P = zeros(ny,nz);
    eig_Vs_GSI_P = zeros(ny,nz);
    eig_Vc_GSI_P = zeros(ny,nz);
    eig_Ws_GSI_P = zeros(ny,nz);
    eig_Wc_GSI_P = zeros(ny,nz);

    eig_Us_GSI = zeros(ny,nz);
    eig_Uc_GSI = zeros(ny,nz);
    eig_Vs_GSI = zeros(ny,nz);
    eig_Vc_GSI = zeros(ny,nz);
    eig_Ws_GSI = zeros(ny,nz);
    eig_Wc_GSI = zeros(ny,nz);

    eig_Us_SFG = zeros(ny,nz);
    eig_Uc_SFG = zeros(ny,nz);
    eig_Vs_SFG = zeros(ny,nz);
    eig_Vc_SFG = zeros(ny,nz);
    eig_Ws_SFG = zeros(ny,nz);
    eig_Wc_SFG = zeros(ny,nz);
    
    C_X = squeeze(CST(:,:,iix,:,:));

    %%
    parfor iiz = 1:Nz_mod
        
        C_XZ = C_X(:,:,:,iiz);

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

        %%%%% Regularization term for eliminating the divergence
        gamma_U = 10;     % Regularization coefficient
        
        Ez = kz(iiz)*eye(ny);
        Ex = diag( -ft(iix)/U )*eye(ny);
        % ss
        ExR = Ex;
        EzR = Ez;
        D3_SS = [ExR , Dy , EzR].*sqrt(d_y);
        IE3 = (eye(3*ny) + gamma_U*D3_SS'*D3_SS);
        inv_GSI_SS = IE3\eye(3*ny);
        % sc
        ExR = Ex;
        EzR = -Ez;
        D3_SC = [ExR , Dy , EzR].*sqrt(d_y);
        IE3 = (eye(3*ny) + gamma_U*D3_SC'*D3_SC);
        inv_GSI_SC = IE3\eye(3*ny);
        % cs
        ExR = -Ex;
        EzR = Ez;
        D3_CS = [ExR , Dy , EzR].*sqrt(d_y);
        IE3 = (eye(3*ny) + gamma_U*D3_CS'*D3_CS);
        inv_GSI_CS = IE3\eye(3*ny);
        % cc
        ExR = -Ex;
        EzR = -Ez;
        D3_CC = [ExR , Dy , EzR].*sqrt(d_y);
        IE3 = (eye(3*ny) + gamma_U*D3_CC'*D3_CC);
        inv_GSI_CC = IE3\eye(3*ny);

        for iiy = 1:Nz_mod

            sy = sin(kz(iiy)*y);
            cy = cos(kz(iiy)*y);
            
            C_XYZ = C_XZ(:,:,iiy);

            %%% The bases for synthesizing the inflow turbulence
            % Part 1
            [eig_U0,eig_SU0] = eig(C_XYZ);

            eig_SU = diag(abs(eig_SU0));    eig_U1 = (eig_U0).*sqrt(eig_SU)';

            [~,sq_eigU] = sort(eig_SU); eig_U2 = eig_U1(:,flipud(sq_eigU));

            eig_Up = [  sy*eig_U2(1,:) ,   cy*eig_U2(1,:)];
            eig_Vp = [  sy*eig_U2(2,:) ,   cy*eig_U2(2,:)];
            eig_Wp = [  sy*eig_U2(3,:) ,   cy*eig_U2(3,:)];
            
            % Part 2
            [eig_U0,eig_SU0] = eig(C_XYZ);

            eig_SU = diag(abs(eig_SU0));    eig_U1 = (eig_U0).*sqrt(eig_SU)';

            [~,sq_eigU] = sort(eig_SU); eig_U2 = eig_U1(:,flipud(sq_eigU));

            eig_Um = [  sy*eig_U2(1,:) ,   cy*eig_U2(1,:)];
            eig_Vm = [  sy*eig_U2(2,:) ,   cy*eig_U2(2,:)];
            eig_Wm = [  sy*eig_U2(3,:) ,   cy*eig_U2(3,:)];
            
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
            phase_init = [ squeeze(phase_random_c(iix,iiz,iiy,:)) ; squeeze(phase_random_s(iix,iiz,iiy,:)) ]' ;
            [phase_opt,J,J_init,hst_1] = f_GSI_1D(eig_Up,eig_Vp,eig_Wp,eig_Um,eig_Vm,eig_Wm,ft(iix),kz(iiz),U,Dy,d_y,phase_init(:));
            
            %%% Synthesize the velocity
            phase_xp = phase_opt(1,:);
            phase_yp = phase_opt(1,:);
            phase_zp = phase_opt(1,:);

            phase_xm = phase_opt(2,:);
            phase_ym = phase_opt(2,:);
            phase_zm = phase_opt(2,:);

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
            Eu_P = (uu_ss.^2 + uu_sc.^2 + uu_cs.^2 + uu_cc.^2)/4;
            Ev_P = (vv_ss.^2 + vv_sc.^2 + vv_cs.^2 + vv_cc.^2)/4;
            Ew_P = (ww_ss.^2 + ww_sc.^2 + ww_cs.^2 + ww_cc.^2)/4;
            S_u_GSI_P(iix,iiy,iiz) = sum(Eu_P*d_y(1));
            S_v_GSI_P(iix,iiy,iiz) = sum(Ev_P*d_y(1));
            S_w_GSI_P(iix,iiy,iiz) = sum(Ew_P*d_y(1));
            
            %%%%% GSI step2: Velocity adjustment %%%%%
            
            % Adjust the eigenmodes to eliminate divergence
            uvw_ss = inv_GSI_CC*([uu_cs(:);vv_ss(:);ww_sc(:)]);
            uvw_sc = inv_GSI_CS*([uu_cc(:);vv_sc(:);ww_ss(:)]);
            uvw_cs = inv_GSI_SC*([uu_ss(:);vv_cs(:);ww_cc(:)]);
            uvw_cc = inv_GSI_SS*([uu_sc(:);vv_cc(:);ww_cs(:)]);

            uu2_cs = uvw_ss(1:ny);vv2_ss = uvw_ss(ny+1:2*ny);ww2_sc = uvw_ss(2*ny+1:end);
            uu2_cc = uvw_sc(1:ny);vv2_sc = uvw_sc(ny+1:2*ny);ww2_ss = uvw_sc(2*ny+1:end);
            uu2_ss = uvw_cs(1:ny);vv2_cs = uvw_cs(ny+1:2*ny);ww2_cc = uvw_cs(2*ny+1:end);
            uu2_sc = uvw_cc(1:ny);vv2_cc = uvw_cc(ny+1:2*ny);ww2_cs = uvw_cc(2*ny+1:end);

            Eu_step2 = sum( (uu2_ss.^2 + uu2_sc.^2 + uu2_cs.^2 + uu2_cc.^2)/4 .* d_y(1));
            Ev_step2 = sum( (vv2_ss.^2 + vv2_sc.^2 + vv2_cs.^2 + vv2_cc.^2)/4 .* d_y(1));
            Ew_step2 = sum( (ww2_ss.^2 + ww2_sc.^2 + ww2_cs.^2 + ww2_cc.^2)/4 .* d_y(1));

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
            S_u_GSI(iix,iiy,iiz) = Eu_step2;
            S_v_GSI(iix,iiy,iiz) = Ev_step2;
            S_w_GSI(iix,iiy,iiz) = Ew_step2;

            %% Reference case: SFG (without divergence-eliminating operations)
            
            %%% Assemble the eigenmodes with random phases
            phase_xp = permute(phase_random_c(iix,iiz,iiy,:),[1 4 2 3]);
            phase_yp = permute(phase_random_c(iix,iiz,iiy,:),[1 4 2 3]);
            phase_zp = permute(phase_random_c(iix,iiz,iiy,:),[1 4 2 3]);

            phase_xm = permute(phase_random_s(iix,iiz,iiy,:),[1 4 2 3]);
            phase_ym = permute(phase_random_s(iix,iiz,iiy,:),[1 4 2 3]);
            phase_zm = permute(phase_random_s(iix,iiz,iiy,:),[1 4 2 3]);

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
            S_u_SFG(iix,iiy,iiz) = sum( (uuR_ss.^2 + uuR_sc.^2 + uuR_cs.^2 + uuR_cc.^2)/4 .* d_y(1));
            S_v_SFG(iix,iiy,iiz) = sum( (vvR_ss.^2 + vvR_sc.^2 + vvR_cs.^2 + vvR_cc.^2)/4 .* d_y(1));
            S_w_SFG(iix,iiy,iiz) = sum( (wwR_ss.^2 + wwR_sc.^2 + wwR_cs.^2 + wwR_cc.^2)/4 .* d_y(1));

        end

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

    %% Output the eigenmodes for DNS computation
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

%%% Energy spectra
E_u_GSI_P = squeeze(sum(S_u_GSI_P(:,1:Nz_mod,1:Nz_mod),[1]));
E_v_GSI_P = squeeze(sum(S_v_GSI_P(:,1:Nz_mod,1:Nz_mod),[1]));
E_w_GSI_P = squeeze(sum(S_w_GSI_P(:,1:Nz_mod,1:Nz_mod),[1]));

E_u_GSI = squeeze(sum(S_u_GSI(:,1:Nz_mod,1:Nz_mod),[1]));
E_v_GSI = squeeze(sum(S_v_GSI(:,1:Nz_mod,1:Nz_mod),[1]));
E_w_GSI = squeeze(sum(S_w_GSI(:,1:Nz_mod,1:Nz_mod),[1]));

E_u_SFG = squeeze(sum(S_u_SFG(:,1:Nz_mod,1:Nz_mod),[1]));
E_v_SFG = squeeze(sum(S_v_SFG(:,1:Nz_mod,1:Nz_mod),[1]));
E_w_SFG = squeeze(sum(S_w_SFG(:,1:Nz_mod,1:Nz_mod),[1]));

[kxx,kzz] = meshgrid(kz(1:Nz_mod),kz(1:Nz_mod));

figure('Position', [500, 500, 800, 600]);

data = {E_u_GSI, E_u_GSI_P, E_u_SFG, E_v_GSI_P, E_v_GSI, E_v_SFG};
titles = {'Euu (GSI-P\&V)', 'Euu (GSI-P)', 'Euu (SFG)', ...
          'Evv (GSI-P\&V)', 'Evv (GSI-P)', 'Evv (SFG)'};

for i = 1:6
    subplot(2, 3, i);
    pcolor(kxx, kzz, data{i} .* kxx .* kzz);
    shading interp;
    set(gca, 'xscale', 'log', 'yscale', 'log');
    title(titles{i}, 'Interpreter', 'latex');
    xlabel('$k_x$', 'Interpreter', 'latex');
    ylabel('$k_z$', 'Interpreter', 'latex');
end

savefig('Energy_spectra_Cubic.fig')

%%% Instantaneous velocity field
[zz,yy] = meshgrid(z,y);

figure('Position',[200,200,1000,600]);

var_data = {uu_GSI(:,:,1),   vv_GSI(:,:,1),   ww_GSI(:,:,1),   duu_GSI+dvv_GSI+dww_GSI, ...
            uu_GSI_P(:,:,1), vv_GSI_P(:,:,1), ww_GSI_P(:,:,1), duu_GSI_P+dvv_GSI_P+dww_GSI_P, ...
            uu_SFG(:,:,1),   vv_SFG(:,:,1),   ww_SFG(:,:,1),   duu_SFG+dvv_SFG+dww_SFG};
        
titles = {'U (GSI-P\&V)', 'V (GSI-P\&V)', 'W (GSI-P\&V)', 'DIV (GSI-P\&V)', ...
          'U (GSI-P)', 'V (GSI-P)', 'W (GSI-P)', 'DIV (GSI-P)', ...
          'U (SFG)', 'V (SFG)', 'W (SFG)', 'DIV (SFG)'};

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
    pcolor(zz,yy,var_data{i});
    shading interp;
    daspect([1 1 1]);
    colorbar;
    caxis(caxis_limits{i});
    title(titles{i}, 'Interpreter', 'latex');
    xlabel('$z$', 'Interpreter', 'latex');
    ylabel('$y$', 'Interpreter', 'latex');
    colormap(jet)
end

savefig('Instantaneous_field_Cubic.fig')

