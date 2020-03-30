% 'test_PnP_with_FFDNet.m' tests Plug-and-Play framework using deep denosing priors (FFDNet) 
% for video reconstruction in 'coded aperture compressive temporal imaging (CACTI)'

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

datasetdir = './dataset';                % dataset dictionary
para.dataname = 'meas_waterBalloon_cr_10';    % selected 2D measurement
para.cr = str2double(para.dataname(end-1:end));    % compression ratio of the selected measurement, i.e., number of video frames to be recovered from each single 2D measurement
para.numRec = 1;                              % number of measurement frames to be reconstructed

datapath = sprintf('%s/%s.mat',datasetdir,para.dataname);  % path of the selected 2D measurement



%% [1] load dataset

load(datapath);               % load measurement
load('./dataset/mask.mat');   % load mask

meas = meas(:,:,1:para.numRec);
meas = para.cr/2*meas./max(meas(:));
mask = double(mask(:,:,1:para.cr));

%% [2] reconstruction

[row, col, ~] = size(meas);
Phi = mask;
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

para.sigma   = [70 60 50]./255; % noise deviation (to be estimated and adapted)
para.vrange   = 1; % range of the signal
para.maxiter = [50 50 50];
para.net = vl_simplenn_tidy(net);

recon = zeros([row, col, para.cr*para.numRec]);
for i_meas = 1:para.numRec
    y = meas(:,:,i_meas);
    recon(:,:,(i_meas-1)*para.cr+(1:para.cr))    =   TV_ADMM_CACTI_FFD_real_GPU( y, para, A,At);
%     recon(:,:,(i_meas-1)*para.cr+(1:para.cr))    =   TV_ADMM_CACTI_FFD_real_GPU_gaptv_ini( y, para, A,At);
end


%% [3] show results in figure

% [3.0] rotate and crop
recon_rotate = zeros(725,725,para.numRec*para.cr);

for np=1:para.numRec*para.cr
    recon_rotate(:,:,np) = imrotate(recon(:,:,np),-135);
end

recon_rotate = recon_rotate(182:182+363,182:182+363,:);
recon_rotate = recon_rotate/max(recon_rotate(:));

% [3.1] show results in figure
figure;
for i=1:para.numRec*para.cr
    imshow(recon_rotate(:,:,i));
    pause(0.2);
end

