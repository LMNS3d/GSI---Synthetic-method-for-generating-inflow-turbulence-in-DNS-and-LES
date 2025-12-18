% Reference paper: Ying, Li & Fu, Generalized synthetic inflow generation %
% method for divergence-free inhomogeneous turbulence, JCP, 2025          %
% Authors: Anjia Ying, Zhigang Li, Lin Fu                                 %

%%%%% Example program: Generate temporal series at the inlet boundary of 
%%%%% TBL using eigenmodes constructed by GSI
%%%%% Prerequisite: Eigenmodes stored in "Mode" obtained by run
%%%%% "main_GSI.m".

%%%%% Note: This program is for example only. In the CFD simulation cases
%%%%% as described in the paper, the velocity signals are generated in situ
%%%%% during the simulation.

clear;close all;
addpath('../Lib_GSI');

dir_ModeData = './Mode/';   % Direction where the divergence-free eigenmodes are stored
dir_InletData = './Inlet/'; % Direction where the generated inflow data are stored
mkdir(dir_InletData)

%% Mean quantities of velocities
load('./Mean_U_V.mat','Umean','Vmean','UU','VV','WW','UV')

Uinf = Umean(end);

%% Grid settings (consistent with the CFD computation, as well as the main_TRL.m)
Nt = 400;
Ny = 151;
Nz = 150;

Lx = 200;
Lz = 15;

t = (1:Nt)/Nt*Lx / Uinf;
z = (1:Nz)/Nz*Lz;

y = load('./yp.dat');

%% Physical time step
Nt_CFD = 10000;
d_tn = 0.002;
tn = (1:Nt_CFD)*d_tn;

%% Wavenumbers and frequencies
Nkz = Nz;
Nft = Nt;

d_t = t(2)-t(1);
d_z = z(2)-z(1);

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

Nzs = Nz/2-1;
Nts = Nt/2-1;

%%
eig_Us = zeros(1,Ny*Nz,Nts);
eig_Uc = zeros(1,Ny*Nz,Nts);
eig_Vs = zeros(1,Ny*Nz,Nts);
eig_Vc = zeros(1,Ny*Nz,Nts);
eig_Ws = zeros(1,Ny*Nz,Nts);
eig_Wc = zeros(1,Ny*Nz,Nts);

for iix = 1:Nts

    fid=fopen([dir_ModeData,'uuS_GSI_iX-',num2str(iix),'.bin']);
    aaa = reshape(fread(fid,'float32'),[Ny,Nz]);
    eig_Us(1,:,iix) = aaa(:);
    fclose(fid);

    fid=fopen([dir_ModeData,'uuC_GSI_iX-',num2str(iix),'.bin']);
    aaa = reshape(fread(fid,'float32'),[Ny,Nz]);
    eig_Uc(1,:,iix) = aaa(:);
    fclose(fid);

    fid=fopen([dir_ModeData,'vvS_GSI_iX-',num2str(iix),'.bin']);
    aaa = reshape(fread(fid,'float32'),[Ny,Nz]);
    eig_Vs(1,:,iix) = aaa(:);
    fclose(fid);

    fid=fopen([dir_ModeData,'vvC_GSI_iX-',num2str(iix),'.bin']);
    aaa = reshape(fread(fid,'float32'),[Ny,Nz]);
    eig_Vc(1,:,iix) = aaa(:);
    fclose(fid);

    fid=fopen([dir_ModeData,'wwS_GSI_iX-',num2str(iix),'.bin']);
    aaa = reshape(fread(fid,'float32'),[Ny,Nz]);
    eig_Ws(1,:,iix) = aaa(:);
    fclose(fid);

    fid=fopen([dir_ModeData,'wwC_GSI_iX-',num2str(iix),'.bin']);
    aaa = reshape(fread(fid,'float32'),[Ny,Nz]);
    eig_Wc(1,:,iix) = aaa(:);
    fclose(fid);
    
end

%% Synthesize and output the inflow velocity with GSI

% Synthesize the velocity signal
ftn = permute(ft(1:Nts),[3 2 1]);

sx = sin(ftn.*tn');
cx = cos(ftn.*tn');

uu_GSI = reshape(sum(( sx.*eig_Us + cx.*eig_Uc ) , 3),[Nt_CFD,Ny,Nz])+Umean';
vv_GSI = reshape(sum(( sx.*eig_Vs + cx.*eig_Vc ) , 3),[Nt_CFD,Ny,Nz])+Vmean';
ww_GSI = reshape(sum(( sx.*eig_Ws + cx.*eig_Wc ) , 3),[Nt_CFD,Ny,Nz]);

disp('Synthesizing velocity signals completed.')

% Write the velocity signal to prescribed direction
uvw_GSI = [uu_GSI(:);vv_GSI(:);ww_GSI(:)];
filename = [dir_InletData,'GSI_inflow'];
fileID = fopen(filename, 'wb');
fwrite(fileID, uvw_GSI, 'double');
fclose(fileID);

disp('Writing velocity signals completed.')

%% Test the results
%%% Instantaneous divergence field
Dy = f_Diff(y,4,'no-slip wall','free');
Dz = f_Diff(z',4,'periodic','periodic');

div_GSI_U = -squeeze(uu_GSI(3,:,:) - uu_GSI(1,:,:))/2/d_tn./Umean;
div_GSI_V = Dy*squeeze(vv_GSI(2,:,:));
div_GSI_W = squeeze(ww_GSI(2,:,:))*Dz';

[zz0,yy0] = meshgrid(z,y);
figure('Position', [500, 500, 1000, 300]);

data = {div_GSI_U + div_GSI_V + div_GSI_W, ...
        -div_GSI_U + div_GSI_V + div_GSI_W};
titles = {'Divergence (GSI)', 'Divergence (reference)'};

for i = 1:2
    subplot(1, 2, i);
    pcolor(zz0, yy0, data{i});
    shading flat;
    colorbar;
    caxis([-1 1]);
    ylim([0 15])
    title(titles{i}, 'Interpreter', 'latex');
    xlabel('$z$', 'Interpreter', 'latex');
    ylabel('$y$', 'Interpreter', 'latex');
end

savefig('TimeSeries_divergence_TBL.fig')

%%% Instantaneous velocity field
[xx1,yy1] = meshgrid(tn*Uinf,y);
[xx2,zz2] = meshgrid(tn*Uinf,z);

figure('Position',[200,200,2000,600]);

data = {squeeze(uu_GSI(:,:,10))', ...
        squeeze(vv_GSI(:,:,10))', ...
        squeeze(uu_GSI(:,10,:))', ...
        squeeze(uu_GSI(:,50,:))'};
xy = {{xx1, yy1}, {xx1, yy1}, {xx2, zz2}, {xx2, zz2}};
titles = {'$u$ (at $t$-$y$ plane)', ...
          '$v$ (at $t$-$y$ plane)', ...
          ['$u$ (at $t$-$z$ plane with $y = ', num2str(y(10)), '$)'], ...
          ['$u$ (at $t$-$z$ plane with $y = ', num2str(y(50)), '$)']};
labels = {{'$t \cdot U$', '$y$'}, ...
          {'$t \cdot U$', '$y$'}, ...
          {'$t \cdot U$', '$z$'}, ...
          {'$t \cdot U$', '$z$'}};
subs = [2,2,1; 2,2,2; 2,2,3; 2,2,4];

for i = 1:4
    subplot(subs(i,1),subs(i,2),subs(i,3));
    pcolor(xy{i}{1}, xy{i}{2}, data{i});
    shading flat;
    colorbar;
    colormap(jet);
    title(titles{i}, 'Interpreter', 'latex');
    xlabel(labels{i}{1}, 'Interpreter', 'latex');
    ylabel(labels{i}{2}, 'Interpreter', 'latex');
end

savefig('TimeSeries_velocity_TBL.fig')

