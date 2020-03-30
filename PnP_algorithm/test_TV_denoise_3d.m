% clear all
% close all
% clc

% Test 3D TV denosing

load('Traffic.mat');
X0 = X(:,:,1:4);

Noise = 0.1* randn(size(X0));
Xnoise = X0 + Noise;

figure; 
subplot(2,2,1); imshow(X0(:,:,1));
subplot(2,2,2); imshow(Xnoise(:,:,1));
psnr(Xnoise, X0)

lambda = 0.15;
iter = 10;

%%
Xde_1 = TV_denoising_clip(Xnoise,lambda,iter);
psnr(Xde_1,X0)
rho = 1e-2;
iter_out = 20; iter_in = 10;
Xde_2 = TV_denoising_ClipB(Xnoise,lambda,rho,iter_out, iter_in);
%Xde_2 = Xde_2./max(Xde_2(:));
psnr(Xde_2,X0)
%figure; 
subplot(2,2,3); imshow(Xde_1(:,:,1));
subplot(2,2,4); imshow(Xde_2(:,:,1));


Xde_3 = tvdenoise_cham_ATV2D(Xnoise,20,iter*2);
psnr(Xde_3, X0)

Xde_4 = fgp_denoise_ATV2D(Xnoise,0.05,40);
psnr(Xde_4, X0)

Xde_5 = tvdenoise_cham_ITV2D(Xnoise,20,iter*2);
psnr(Xde_5,X0)

Xde_6 = tvdenoise_cham_ITV3D(Xnoise,20,iter*2);
psnr(Xde_6,X0)

Xde_7 = fgp_denoise_ITV2D(Xnoise,0.05,40);
psnr(Xde_7, X0)

Xde_8 = fgp_denoise_ITV3D(Xnoise,0.1,40); 
psnr(Xde_8, X0)