% 'test_DeSCI.m' tests 'decompress snapshot compressive imaging (DeSCI)' algorithm for video reconstruction in
% 'coded aperture compressive temporal imaging (CACTI)'

% Reference
%   [1] M. Qiao, Z. Meng, J. Ma, X. Yuan, Deep learning for video compressive
%       sensing, APL Photonics 5, 030801 (2020).
%   [2] Y. Liu, X. Yuan, J. Suo, D.J. Brady, and Q. Dai, Rank Minimization 
%       for Snapshot Compressive Imaging, IEEE Trans. Pattern Anal. Mach. 
%       Intell. (TPAMI), DOI:10.1109/TPAMI.2018.2873587, 2018.


% Contact
%   Xin Yuan, Bell Labs, xyuan@bell-labs.com
%   Mu Qiao, New Jersey Institute of Technology, muqiao@njit.edu
%   Update Mar 13, 2020.


%% [0] environment configuration
clear; 
clc;
close all

addpath(genpath('./DeSCI_algorithm')); % algorithms
datasetdir = './dataset';                % dataset dictionary

para.dataname = 'meas_waterBalloon_cr_10';    % selected 2D measurement
para.cr = str2double(para.dataname(end-1:end));    % compression ratio of the selected measurement, i.e., number of video frames to be recovered from each single 2D measurement
para.numRec = 1;                              % number of measurement frames to be reconstructed

datapath = sprintf('%s/%s.mat',datasetdir,para.dataname);  % path of the selected 2D measurement

%% [1] load dataset

load(datapath);               % load measurement
load('./dataset/mask.mat');   % load mask

meas = meas(:,:,1:para.numRec);
meas = 1850*meas./max(meas(:));
mask  = double(mask(:,:,1:para.cr));

[nrow,ncol,~] = size(meas);

%% [2] ADMM-WNNM-TV
para.nframe =   para.numRec; 
para.MAXB   = 255;
MAXB = para.MAXB;

para.Mfunc  = @(z) A_xy(z,mask);
para.Mtfunc = @(z) At_xy_nonorm(z,mask);
para.Phisum = sum(mask.^2,3);
para.Phisum(para.Phisum==0) = 1;

para.flag_iqa = false; % disable image quality assessments in iterations
para.acc = 1; % enable acceleration
para.flag_iqa = false; % disable image quality assessments in iterations
para.projmeth = 'admm_res'; % projection method
%  (GAP for noiseless or ADMM for noisy)
para.gamma  =    1; % regularization factor for noise suppression
para.denoiser = 'wnnm'; % WNNM denoising
para.wnnm_int_fwise = true; % enable GAP-WNNM integrated (with frame-wise denoising)
para.blockmatch_period = 20; % period of block matching
para.sigma   = [50 45 40]/MAXB; % noise deviation (to be estimated and adapted)
para.vrange   = 1; % range of the signal
para.maxiter = [50 50 50]; % first trial
para.patchsize = 24; % patch size
para.iternum = 1; % iteration number in WNNM
para.enparfor = true; % enable parfor
para.numworkers = 12;
if para.enparfor % if parfor is enabled, start parpool in advance
    mycluster = parcluster('local');
    delete(gcp('nocreate')); % delete current parpool
    ord = 0;
    while para.cr/2^ord > mycluster.NumWorkers
        ord = ord+1;
    end
    poolobj = parpool(mycluster,min(max(floor(para.cr/2^ord),1),para.numworkers));
end

[vadmmwnnmtv,~,~,tadmmwnnm] = ...
    admmdenoise_cacti(mask,meas,[],[],para);

%% [3] show results in figure

% [3.0] rotate and crop
recon = vadmmwnnmtv;
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

