% generate PD, T2s, T1 phantoms from chi phantom at 3 Tesla

% chi phanotm originally provided by:
% Milovic C, et al. FAst Nonlinear Susceptibility Inversion (FANSI Toolbox for QSM)
%, 34th ESMRMB Annual Scientific Meeting. BARCELONA

%% chnage resolution to 0.868x0.868x0.85 mm
res = [256*1 256*1 98*1.5]./[295 295 173];

%% laod masks
load brain_masks


 
%% PD map
%        1.CD   2.CSF   3.GM    4.GP      5.PT    6.SN    7.TH     8.WM
PD_w =  [.802,  1.00,   .807   ,.744 ,   .797 ,  0.75,   .756 ,    .679 ];

PD_map =  PD_w(1)*cd_m + PD_w(2)*csf_m +PD_w(3)*gm_m + PD_w(4)*gp_m +PD_w(5)*pt_m ...
    + PD_w(6)*sn_m +PD_w(7)*th_m + PD_w(8)*wm_m ;


% nii = make_nii(PD_map,res);
%  save_nii(nii,'PD_map.nii');


 %% Susceptibility map
%        1.CD   2.CSF   3.GM    4.GP      5.PT    6.SN    7.TH     8.WM
sus_w = [0.03,   0.0,  -0.01,   0.10,     0.04,    0.09 ,  0.0,   -0.04];

sus_map = sus_w(1)*cd_m + sus_w(2)*csf_m +sus_w(3)*gm_m + sus_w(4)*gp_m +sus_w(5)*pt_m ...
    + sus_w(6)*sn_m +sus_w(7)*th_m + sus_w(8)*wm_m ;

%  nii = make_nii(sus_map,res);
%  save_nii(nii,'sus_map.nii');
 
 %% T2* map
 %        1.CD   2.CSF   3.GM    4.GP      5.PT    6.SN    7.TH     8.WM
 T2s_w = [54.8,   2000,  58,     28.8,     53.6,   42,    76.8,      55];

 
 T2s_map =  T2s_w(1)*cd_m + T2s_w(2)*csf_m +T2s_w(3)*gm_m + T2s_w(4)*gp_m +T2s_w(5)*pt_m ...
    + T2s_w(6)*sn_m +T2s_w(7)*th_m + T2s_w(8)*wm_m ;


% nii = make_nii(T2s_map,res);
%  save_nii(nii,'T2s_map.nii');

 
 %% T1 map
 %        1.CD   2.CSF   3.GM    4.GP      5.PT    6.SN    7.TH     8.WM
 T1_w = [1280,   4500,   1390,   950,      1150,   925,    1150,    910];  

 
  T1_map =   T1_w(1)*cd_m + T1_w(2)*csf_m +T1_w(3)*gm_m + T1_w(4)*gp_m +T1_w(5)*pt_m ...
    + T1_w(6)*sn_m +T1_w(7)*th_m + T1_w(8)*wm_m ;


% nii = make_nii(T1_map,res);
%  save_nii(nii,'T1_map.nii');
 
 
 
 %% reduce FOV in PE dir to 208 lines
 
  PD_map = PD_map(44:end-44,:,:);
  sus_map = sus_map(44:end-44,:,:);
  T2s_map = T2s_map(44:end-44,:,:);
  T1_map = T1_map(44:end-44,:,:);
  
  cd_m = cd_m(44:end-44,:,:);
  pt_m = pt_m(44:end-44,:,:);
  th_m = th_m(44:end-44,:,:);
  gp_m = gp_m(44:end-44,:,:); 
  wm_m = wm_m(44:end-44,:,:);
  csf_m = csf_m(44:end-44,:,:);
  sn_m = sn_m(44:end-44,:,:);    
  gm_m = gm_m(44:end-44,:,:);   

 
  mask = cd_m + pt_m + th_m + gp_m + wm_m + csf_m + sn_m + gm_m;
  mask(mask~=0) = 1;
%% show maps
sel = 85;

scrsz = get(groot,'ScreenSize');
figure('Position',[scrsz(3)/2 scrsz(4)/3 scrsz(3)/3.5 scrsz(4)/2])
subplot(2,2,1),imagesc(rot90(PD_map(:,:,sel),1)),title('PD map '),caxis([min(PD_map(:)) max(PD_map(:))]),colorbar,colormap gray
subplot(2,2,2),imagesc(rot90(sus_map(:,:,sel),1)),title('sus map (ppm) '),caxis([-0.10 0.10]),colorbar,colormap gray
subplot(2,2,3),imagesc(rot90(T1_map(:,:,sel),1)),title('T1 map (ms) '),caxis([min(T1_map(:)) max(T1_map(:))]),colorbar,colormap gray
subplot(2,2,4),imagesc(rot90(T2s_map(:,:,sel),1)),title('T2* map (ms) '),caxis([min(T2s_map(:)) 100]),colorbar,colormap gray   
     
%% save
clear sel scrsz
save('brain_maps_208.mat')