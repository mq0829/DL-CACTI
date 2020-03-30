%% para
% 
cr = 10;    % overall cr

filename = ['waterBalloon_cr_' num2str(cr)];

numRec = 1;                % number of coded frames to be reconstructed

sigma = [70];

maxiter = 5*ones(1,length(sigma));
maxiter = [50];

%%

load(['meas_',filename,'.mat']);
load('mask.mat');
load(['recon_gaptv_',filename,'.mat']);

% meas = meas(:,:,1);
meas = meas./max(meas(:))*cr/2;
mask = mask(:,:,1:cr);
global recon_ini;   recon_ini = recon;   delete recon;
recon_ini = recon_ini(:,:,1:cr);

%%

[row, col, numMeas] = size(meas);
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
para.iter =150;
para.acc= 1;
para.mu = 0.25;

% para.maxiter = [10 10 10 10]; % first trial

TV_m_all = {'ATV_ClipA', 'ATV_ClipB','ATV_cham','ATV_FGP','ITV2D_cham','ITV2D_FGP','ITV3D_cham','ITV3D_FGP'};

% addpath(genpath('D:\learning_Code\FFDNet-master\FFDNet-master'));

para.eta=  0.5;

para.TVm = TV_m_all{1};

%% now we use FDDnet
para.TV_iter = 1;
para.sigma   = sigma./255; % noise deviation (to be estimated and adapted)
para.vrange   = 1; % range of the signal
para.maxiter = maxiter;

load(fullfile('.\FFD_net_model','FFDNet_gray.mat'));
addpath(genpath('./funs'));
addpath(fullfile('./FFDNet-master/utilities'));

para.net = vl_simplenn_tidy(net);

tic
vffenet = zeros([row, col, cr*numMeas]);
for i_meas = 1:1
    y = meas(:,:,i_meas);
%     vffenet(:,:,(i_meas-1)*cr+(1:cr))    =   TV_ADMM_CACTI_FFD_real_GPU( y, para, A,At);
    vffenet(:,:,(i_meas-1)*cr+(1:cr))    =   TV_ADMM_CACTI_FFD_real_GPU_gaptv_ini( y, para, A,At);
end
toc

%% [3] show and save results

recon = vffenet/max(vffenet(:)); % normalization
recon(recon<0)=0;

% show results in figure
for i=5
    figure;imshow(recon(:,:,i));
    title([num2str(para.sigma*255) ',  ' num2str(para.maxiter)]);
  
%     figure;imshow(recon_ini(:,:,i));
end


% save results to .mat file
% save(['recon_ffdnet_',filename,'.mat'],'recon');


