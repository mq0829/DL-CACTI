% load and show various TV
clear all
close all
clc

for nfi = 1:5
switch nfi
    case 1
      fname = 'kobe';
    case 2
      fname = 'traffic';
    case 3 
      fname = 'runner';   
    case 4
      fname = 'drop';
    case 5
      fname = 'crash';
end
load([fname '_all4_variousTV_20190930.mat'])

for nn=1:8
    psnr_fista(nn,nfi) = psnr(theta_fista{nn}, orig./255);
    psnr_twist(nn,nfi) = psnr(theta_twist{nn}, orig./255);
    psnr_gap(nn,nfi) = psnr(theta_gap{nn}, orig./255);
    psnr_admm(nn,nfi) = psnr(theta_admm{nn}, orig./255);
end

end

psnr_fista_mean = mean(psnr_fista, 2);
psnr_twist_mean = mean(psnr_twist, 2);
psnr_gap_mean = mean(psnr_gap, 2);
psnr_admm_mean = mean(psnr_admm, 2);

psnr_fista_mean4 = mean(psnr_fista(:,1:4), 2);
psnr_twist_mean4 = mean(psnr_twist(:,1:4), 2);
psnr_gap_mean4 = mean(psnr_gap(:,1:4), 2);
psnr_admm_mean4 = mean(psnr_admm(:,1:4), 2);



  TV_m_all = {'ATV-ClipA', 'ATV-ClipB','ATV-Cham','ATV-FGP','ITV2D-Cham','ITV2D-FGP','ITV3D-cham','ITV3D-FGP'};
  aglo = 'GAP';
  
  tf = 10;
  figure;
  subplot(2,5,1);  imshow(orig(:,:,tf)./255);  title('Truth');
  for np = 1:8
      if(np<5)
      subplot(2,5,np+1);
      else
         subplot(2,5,np+2);
      end
      imshow(theta_gap{np}(:,:,tf));  title([aglo '-'  TV_m_all{np} ', ' num2str(psnr(theta_gap{np}(:,:,tf),orig(:,:,tf)./255)) 'dB']);
  end

tf = 10; 
figure; subplot(2,2,1);  imshow(orig(:,:,tf)./255);  title('Truth');
subplot(2,2,3); imshow(theta_gap{5}(:,:,tf)); title('GAP-ITV2D-Cham');
subplot(2,2,2); imshow(theta_gap{1}(:,:,tf)); title('GAP-ATV-ClipA');
subplot(2,2,4); imshow(theta_gap{8}(:,:,tf));  title('GAP-ITV3D-FGP');

for tf = 1:40
    psnr_atv_clip(tf) = psnr(squeeze(theta_gap{1}(:,:,tf)), squeeze(orig(:,:,tf)./255));
    psnr_itv2d_cham(tf) = psnr(squeeze(theta_gap{5}(:,:,tf)), squeeze(orig(:,:,tf)./255));
    psnr_itv3d_fgp(tf) = psnr(squeeze(theta_gap{8}(:,:,tf)), squeeze(orig(:,:,tf)./255));
end
figure; plot(psnr_atv_clip,'k-'); 
hold on; plot(psnr_itv2d_cham,'r-');
hold on; plot(psnr_itv3d_fgp,'g-')

