% load and read data for Mu Qiao's CACTI using DMD
% date: 2019-10-10

filename = 'labs_Real_20191005';
load([filename '.mat']);

addpath(genpath('./funs'));


[row, col, ColT] = size(mask_bayer);
y = meas_bayer./max(meas_bayer(:))*ColT/2;
Phi = mask_bayer;

para.row = row;
para.col = col;

para.TVweight = 0.1;

A = @(z) A_xy(z, Phi);
%At = @(z) At_xy(z, Phi,Phi_sum);
At = @(z) At_xy_nonorm(z, Phi);

 Phi_sum = sum(Phi.^2,3);
Phi_sum(Phi_sum==0)=1;
 para.lambda = 1;
 para.Phi_sum = Phi_sum;
 para.iter =150;
 para.acc= 1;
 para.mu = 0.25;
 
 

  % para.maxiter = [10 10 10 10]; % first trial
 
  TV_m_all = {'ATV_ClipA', 'ATV_ClipB','ATV_cham','ATV_FGP','ITV2D_cham','ITV2D_FGP','ITV3D_cham','ITV3D_FGP'};
  
 % addpath(genpath('D:\learning_Code\FFDNet-master\FFDNet-master'));

para.eta=  0.5;

para.TVm = TV_m_all{1};
% im    =   TV_ADMM_CACTI_real( y, para, A,At);
% 
% figure;
% for tt=1:ColT
%     subplot(2, ceil(ColT/2),tt); imshow(im(:,:,tt)/max(im(:)));
% end


%% now we use FDDnet
para.TV_iter = 10;
 para.sigma   = [50 40 30 ]./255; % noise deviation (to be estimated and adapted)
  para.vrange   = 1; % range of the signal
  para.maxiter = [ 60 60 60 ]./6;
    
load(fullfile('.\FFD_net_model','FFDNet_gray.mat'));
para.net = vl_simplenn_tidy(net);
 im_ffd    =   TV_ADMM_CACTI_FFD_real( y, para, A,At);
 figure;
for tt=1:ColT
    subplot(2, ceil(ColT/2),tt); imshow(im_ffd(:,:,tt)/max(im_ffd(:)));
end
