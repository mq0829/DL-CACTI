% 'test_DeSCI.m' tests 'Generalized alternating projection based total variation minimization (GAP-TV)' algorithm 
% for video reconstruction in 'coded aperture compressive temporal imaging (CACTI)'

% Reference
%   [1] M. Qiao, Z. Meng, J. Ma, X. Yuan, Deep learning for video compressive
%       sensing, APL Photonics 5, 030801 (2020).
%   [2] Yuan, Xin. "Generalized alternating projection based total variation minimization for compressive sensing." 
%       2016 IEEE International Conference on Image Processing (ICIP). IEEE, 2016.


% Contact
%   Xin Yuan, Bell Labs, xyuan@bell-labs.com
%   Mu Qiao, New Jersey Institute of Technology, muqiao@njit.edu
%   Update Mar 13, 2020.

% The main reconstruction function 'TV4_ADMM_CACTI_adaw' called in line 63 and the subfunctions it invokes are writen fully in Matlab code
% A faster version, 'TV4_ADMM_CACTI_adaw_ap', is aslo available in this package that invokes a C-written denosing core ('tvden.mexw64') 
% and shortens the run-time by ~5X (Matlab2018a or higer version required)

%% [0] environment configuration
clear; 
clc;
close all

addpath(genpath('./GAP_TV_algorithm'));             % algorithm package

datasetdir = './dataset';                           % dataset dictionary
para.dataname = 'meas_waterBalloon_cr_10';          % select a 2D measurement
para.cr = str2double(para.dataname(end-1:end));     % compression ratio of the selected measurement, i.e., number of video frames to be recovered from each single 2D measurement
para.numRec = 1;                                    % number of measurement frames to be reconstructed

datapath = sprintf('%s/%s.mat',datasetdir,para.dataname);  % path of the selected 2D measurement

%% [1] load dataset
load(datapath);               % load measurement
load('./dataset/mask.mat');   % load mask

meas = meas(:,:,1:para.numRec);
meas = 255*meas/max(meas(:));
mask = double(255*mask(:,:,1:para.cr));
mask = mask/max(mask(:));

%% [2] reconstruction

% [2.1] algorithm para
A = @(z) A_xy(z, mask);
At = @(z) At_xy_nonorm(z, mask);

Phi_sum = sum(mask.^2,3);
Phi_sum(Phi_sum==0)=1;
para.lambda = 1;  
para.Phi_sum = Phi_sum;
para.iter =100;     
para.TVweight = 1; 
para.eta = 10; 

% [2.2] run reconstruction
recon = zeros([size(meas,1),size(meas,2),para.cr*size(meas,3)]);
tic;
for i_meas = 1:para.numRec
    meas_single = meas(:,:,i_meas);
    
    vgaptv  =  TV4_ADMM_CACTI_adaw(meas_single, para, A,At);

    recon(:,:,(i_meas-1)*para.cr +(1:para.cr)) = vgaptv;
    
    fprintf( 'ADMM-TV, %dth frame out of %d, time to now: %d s\n', i_meas, para.numRec, round(toc) );
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

