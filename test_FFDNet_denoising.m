% 'test_FFDNet_denoising.m' performs frame-wise denoising on video reconstructions from other algorithms (GAP-TV, DeSCI or PnP) 
% using deep denosing FFDNet 

% Reference
%   [1] M. Qiao, Z. Meng, J. Ma, X. Yuan, Deep learning for video compressive
%       sensing, APL Photonics 5, 030801 (2020).
%   [2] Wang, X., & Chan, S. H. (2017, March). Parameter-free plug-and-play ADMM for image restoration. 
%       In 2017 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP) (pp. 1323-1327). IEEE.


% Contact
%   Xin Yuan, Bell Labs, xyuan@bell-labs.com
%   Mu Qiao, New Jersey Institute of Technology, muqiao@njit.edu
%   Update Mar 13, 2020.

% For environment requirements, refer to 'https://github.com/cszn/FFDNet'

%% [0] environment configuration
clear; 
clc;
close all

addpath(genpath('./PnP_algorithm/funs'));
addpath(fullfile('./PnP_algorithm/FFDNet-master/utilities'));
addpath('./PnP_algorithm\matconvnet-1.0-beta25\matlab\simplenn');

addpath(genpath('.\PnP_algorithm/FFD_net_model'));
load(fullfile('.\PnP_algorithm/FFD_net_model','FFDNet_gray.mat'));
net = vl_simplenn_move(vl_simplenn_tidy(net), 'gpu');

datasetdir = './dataset';                                   % dataset dictionary
para.dataname = 'recon_gaptv_waterBalloon_cr_10';           % select a reconstructed video
para.cr = str2double(para.dataname(end-1:end));             % compression ratio

datapath = sprintf('%s/%s.mat',datasetdir,para.dataname);   % path of the selected file
noise_est = 50/255;                                         % noise level estimate

%% [1] load noisy video reconstruction
load(datapath);
para.numRec = size(recon,3)/para.cr;

%% [2] perform denoising
recon_gaptv = recon; clear recon;
recon_denoise = zeros(size(recon_gaptv));
for ncc = 1:numRec*para.cr
    recon_denoise(:,:,ncc)         =  FFD_Net_DenoiserGPU(single(recon_gaptv(:,:,ncc)), noise_est,net);
end

%% [3] show results in figure

% [3.0] rotate and crop 'recon_denoise'
recon_denoise_rotate = zeros(725,725,para.numRec*para.cr);

for np=1:para.numRec*para.cr
    recon_denoise_rotate(:,:,np) = imrotate(recon_denoise(:,:,np),-135);
end

recon_denoise_rotate = recon_denoise_rotate(182:182+363,182:182+363,:);
recon_denoise_rotate = recon_denoise_rotate/max(recon_denoise_rotate(:));


% [3.1] rotate and crop 'recon_gaptv'
recon_gaptv_rotate = zeros(725,725,para.numRec*para.cr);

for np=1:para.numRec*para.cr
    recon_gaptv_rotate(:,:,np) = imrotate(recon_gaptv(:,:,np),-135);
end

recon_gaptv_rotate = recon_gaptv_rotate(182:182+363,182:182+363,:);
recon_gaptv_rotate = recon_gaptv_rotate/max(recon_gaptv_rotate(:));


% [3.2] compare results
figure;
for i=1:para.numRec*para.cr
    subplot(1,2,1);
    imshow(recon_gaptv_rotate(:,:,i));
    subplot(1,2,2);
    imshow(recon_denoise_rotate(:,:,i));    
    pause(0.2);
end
