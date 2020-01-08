% prdouce MPRAGE mag and phase using brain phantom

function MPRAGE_phase_simulation
% load brain phantom data
load brain_maps_208

B0 = 3;                  % T, Imaging field strength
gyro = 2*pi*42.58*1e6;   % rad/T, gyromagnetic ratio
TR = 1800;               % ms, repitition time inv-to-inv
TE = 2.37;               % ms, echo time
TIg = 900;               % ms, inversion time 

TS = 7.20 ;              % ms, spaceing between readouts 
N = 208;                 % turbo factor, # of readouts in each TR
a = 180;                 % degree, inversion flip-angle
theta = 8;               % degree, readout flip-angle 
tau = N*TS  ;            % total time of all readouts

SNR = 32;                % [46, 21.5, 13.8, 9.5, 2.75, 1.53];

sim_mode = 'user';       % 'user' or 'max TI'

TI_mode = 'middle';      % 'first' or 'middle', assume TI until frist or middle readout

T1 = T1_w;               % ms, T1 longitudinal relaxation time
T2 = T2s_w;              % ms, T2* relaxation time

Nss =  10;               % repeat few cycles (few TRs) to achieve steady response 
%% TI calculation

if strcmp(sim_mode,'user')  % use the specific TI given by user
    
    if strcmp(TI_mode,'middle')
        TI =  TIg- tau/2;
        if TI<0
            error('negative TI !')
        end
    else 
        TI =  TIg;
    end


    TD = TR - TI- tau; % ms, delay time, 

    if TD < 0
        error('delay time is unrelizable!')
    end

    
    
elseif strcmp(sim_mode,'max TI')   % simulate all possible TI
    
    TI_max = TR - tau;  % assuming  TD = 0
    
    if TI_max<0
            error('negative TI !')
    end
    
    
    [T1,TI] = meshgrid(T1, 1:TI_max);
    TD = TR - TI- tau; % ms, delay times 
end
%% bloch simulation 


Mzr = zeros([size(T1) N]); % stores Mz value at TE=0;
M0 = [0; 0; 1];            % inital magnetization along z-axis
Mn = M0;
Mt = Mn(1) + 1i.*Mn(2);    % transverse component of magnetization
Mz = Mn(3);                % longitudinal component of magnetization
l_Mt = Mt;
l_Mz = Mz;

for r =1:Nss  % for Nss TRs
    
  % inversion pulse along x
  Mt = real(l_Mt) + 1i* (imag(l_Mt).*cosd(a) + l_Mz.*sind(a));
  Mz = (-imag(l_Mt).*sind(a) + l_Mz.*cosd(a));
  
  % wait for TI
   Mt = Mt.*exp(-TI./T2); 
   Mz = Mn(3)+(Mz-Mn(3)).*exp(-TI./T1);
   
   l_Mt = Mt;
    l_Mz = Mz;
   
   for n =1:N  % loop for turbo RF pulses
       
       % RF pulse with theta along x
       Mt = real(l_Mt) + 1i* (imag(l_Mt).*cosd(theta) + l_Mz.*sind(theta));
       Mz = (-imag(l_Mt).*sind(theta) + l_Mz.*cosd(theta));
       
       Mzr(:,:,n)= Mz; % store Mz value 

       
       % wait for TS
         Mt = Mt.*exp(-TS./T2); 
         Mz = Mn(3)+(Mz-Mn(3)).*exp(-TS./T1);
        
%          ASSUME Perfect spoiling 
         Mt = 0;
         l_Mt = Mt;
         l_Mz = Mz;
   end
   
   % wait for TD time
   
         Mt = Mt.*exp(-TD./T2); 
         Mz = Mn(3)+(Mz-Mn(3)).*exp(-TD./T1);
   
          % ASSUME Perfect spoiling 
         Mt = 0;
         l_Mt = Mt;
         l_Mz = Mz;
end

%% construct magnitude  & phase
mag = PD_map;
% normalize mag to 1
mag = mag./max(mag(:));

siz = size(mag);

% generate mask
mask = zeros(siz); 
mask(mag~=0) = 1;
mag = mag + 0.01;
%% dipole kernel
D = dipole_kernel(siz, res);
% prepare local field
b1 = real(ifftn(D.*fftn(sus_map)));

% backgurnd field
% b2 = real(ifftn(D.*fftn((1-mask))));

% total field
b = b1;% +b2;

% scale into phase 
phs = b*1e-6 *B0*gyro*TE*1e-3;

% initial complex MR signal
Mz = mag.*exp(1i.*phs) ;

Mzr = permute(Mzr, [2 3 1]);

% include flip-angle effect and t2* relaxation
sig2 = Mzr .* sind(theta) .* repmat(exp(-TE./T2s_w.'),1,size(Mzr,2));
%% construct kspace line by line considering relaxation effect

for i=1:N   % for each k-space line along slice encoding 

    % prepare MR signal in space domain
    sig_r = Mz.*(cd_m.*sig2(1,i) + csf_m.*sig2(2,i)+gm_m.*sig2(3,i) + gp_m.*sig2(4,i) + pt_m.*sig2(5,i) + sn_m.*sig2(6,i)...
         +th_m.*sig2(7,i)+wm_m.*sig2(8,i));  
     
    % transform into k-space
    ks = fftshift(fftn(sig_r));
    
    % store one slice encoding line 
    ks_rr(i,:,:) = ks (i,:,:); % y = slice encoding direction
end

% Pepare noise
noise = (randn(siz)+1i*randn(siz));
n_std = std(noise(:));
noise = noise/n_std;   % normalize noise std to 1

mean(abs(ks_rr(:)))
% add noise in k-space 
ks_r = (ks_rr )+ (noise/SNR*53.5);
% ks_r = elept_filt(ks_r);

% transform into space domain
sign_rcv = ifftn(ifftshift(ks_r));

% store MPRAGE magnitude and phase 
MPRAGE_mag = abs(sign_rcv);
MPRAGE_phs = angle(sign_rcv);
% nii = make_nii(MPRAGE_mag);save_nii(nii,['MPRAGE_mag_',num2str(SNR),'.nii']);
% nii = make_nii(MPRAGE_phs);save_nii(nii,['MPRAGE_phs_',num2str(SNR),'.nii']);


%% measure SNR on magnitdue 
[m_snr,snr_te, men,st] = calc_snr(sign_rcv,gp_m+pt_m+th_m,TE);

m_snr


 %% CNR calc
 CNR_cd_ic = (sus_w(1) - sus_w(8)).*1e-6 *B0*gyro.*snr_te*1/3; % delta chi calculated with ref to WM
 CNR_gp_ic = (sus_w(4) - sus_w(8)).*1e-6 *B0*gyro.*snr_te*1/3;
 CNR_pt_ic = (sus_w(5) - sus_w(8)).*1e-6 *B0*gyro.*snr_te*1/3;
 CNR_th_ic = (sus_w(7) - sus_w(8)).*1e-6 *B0*gyro.*snr_te*1/3; 
 
%% Show one axial image
sel = 85;

scrsz = get(groot,'ScreenSize');
figure('Position',[scrsz(3)/2 scrsz(4)/3 scrsz(3)/3.5 scrsz(4)/2])
subplot(2,2,1),imagesc(rot90(MPRAGE_mag(:,:,sel),1)),title(['MPRAGE magnt, SNRm:',num2str(m_snr,3)]),caxis([0 0.09]),colorbar,colormap gray
subplot(2,2,2),imagesc(rot90(MPRAGE_phs(:,:,sel),1)),title('MPRAGE phase (rad) '),caxis([-.14 .14]),colorbar,colormap gray
subplot(2,2,3:4), plot([1;208],[0;0],':k','LineWidth',1.5) % dashed line for Mz =0
hold on
Mzr2 = Mzr; Mzr2(6:7,:) = []; % no need to show TH and SN

plot(1:N,Mzr2.','LineWidth',2),title(['Mz evolution during readouts, TI: ',num2str(TIg), 'ms, TR: ',num2str(TR),'ms, FA:',num2str(theta),'^{\circ}']),xlim([1 208]),xlabel('Readouts'),ylabel('Mz'),grid on
plot(104,0,'pk','MarkerSize',9,'MarkerFaceColor','k') % k-space center, middle readout
ax = gca;
ax.XTick = [4:50:208];

legend('Mz = 0','CD','CSF','GM','GP','PT & TH','WM','Middle Readout')

 %% save
 clear  pt_m gp_m th_m cd_m sn_m gm_m wm_m csf_m mask T2s_map T1_map PD_map sus_map Mzr2
     
nam =['MPRAGE_TE',num2str(TE),'_SNRm_',num2str(m_snr),'.mat'];

if exist(nam,'file') ==2
    nam = ['MPRAGE_TE',num2str(TE),'_SNRm_',num2str(m_snr),'_2.mat'];
end
 save(nam)
end

function [snr, snr_te, m,st] = calc_snr(sign_rcv,msk,te)
 % MEAN
 m = abs(sign_rcv).*msk; 
 m = sum(m(:))./sum(msk(:));
 % STD in
 st = ((abs(sign_rcv) - m).*msk).^2;
 st = sqrt(sum(st(:))./sum(msk(:)));
 
 % air mask
 msk_noise = zeros(size(msk));
 msk_noise(20:170,1:20,20:150) = 1;
 
 % noise STD
 m_N = abs(sign_rcv).*msk_noise; 
 m_N = sum(m_N(:))./sum(msk_noise(:));
 st_N = ((abs(sign_rcv) - m_N).*msk_noise).^2;
 st_N = sqrt(sum(st_N(:))./sum(msk_noise(:)));
 
 % SNR
 snr = m/st_N/1.526; % correction factor  = sqrt(2/(4-pi))
 snr_te = snr.*te*1e-3;
end

function y = elept_filt(x)
N = size(x);
[ky,kx,kz] = meshgrid(-N(2)/2:N(2)/2-1, -N(1)/2:N(1)/2-1, -N(3)/2:N(3)/2-1);
a = N(1)/2;
b = N(2)/2;
c = N(3)/2;

mask = zeros(N);
mask((kx/a).^2 + (ky/b).^2 + (kz/c).^2 <=1) = 1;
y = x.*mask;
end

function [ kernel ] = dipole_kernel( N, spatial_res )
% input:
% N - array size
% spatial_res - 
% continuous kernel proposed by Salomir, et al. 2003.
% output:
% kernel - dipole kernel in the frequency space
[ky,kx,kz] = meshgrid(-N(2)/2:N(2)/2-1, -N(1)/2:N(1)/2-1, -N(3)/2:N(3)/2-1);
kx = (kx / max(abs(kx(:)))) / spatial_res(1);
ky = (ky / max(abs(ky(:)))) / spatial_res(2);
kz = (kz / max(abs(kz(:)))) / spatial_res(3);

k2 = kx.^2 + ky.^2 + kz.^2;

R_tot = eye(3);

kernel = fftshift( 1/3 - (kx * R_tot(3,1) + ky * R_tot(3,2) + kz * R_tot(3,3)).^2 ./ (k2 + eps) );   
end
