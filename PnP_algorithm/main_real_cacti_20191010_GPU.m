object = {'hand_cr_10','hand_cr_20','hand_cr_30','hand_cr_40','hand_cr_50',...
        'duomino_cr_10','duomino_cr_20','duomino_cr_30','duomino_cr_40','duomino_cr_50',...
        'pendulumBall_cr_10','pendulumBall_cr_20','pendulumBall_cr_30','pendulumBall_cr_40','pendulumBall_cr_50',...
        'pingpang_cr_10','pingpang_cr_20','pingpang_cr_30','pingpang_cr_40','pingpang_cr_50',...
        'waterBalloon_cr_10','waterBalloon_cr_20','waterBalloon_cr_30','waterBalloon_cr_40','waterBalloon_cr_50'};
crs = [10 20 30 40 50 10 20 30 40 50 10 20 30 40 50 10 20 30 40 50 10 20 30 40 50];

for cc = 1:25
%% para

filename = object{cc};
cr = crs(cc);


frame_rate = 5;            % frame rate of the saved video
time_interval = 2*10/cr;        % ms


%%

load(['meas_',filename,'.mat']);
load('mask.mat');

addpath(genpath('./funs'));
addpath(fullfile('./FFDNet-maser/utilities'));

[row, col, numMeas] = size(meas);
Phi = mask(:,:,1:cr);
meas = meas./max(meas(:))*cr/2;

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
para.sigma   = [50 40 30 ]./255; % noise deviation (to be estimated and adapted)
para.vrange   = 1; % range of the signal
para.maxiter = [ 60 60 60 ]./6;

load(fullfile('.\FFD_net_model','FFDNet_gray.mat'));
para.net = vl_simplenn_tidy(net);

tic
im_ffd = zeros([size(meas),cr]);
for i_meas = 1:numMeas
    y = meas(:,:,i_meas);
    im_ffd(:,:,i_meas,:)    =   TV_ADMM_CACTI_FFD_real_GPU( y, para, A,At);
end
toc

%% [3] show and save results

im_ffd_scale = im_ffd/max(im_ffd(:)); % normalization
im_ffd_scale(im_ffd_scale<0)=0;

% show results in figure
for i=1
    figure;imshow(im_ffd_scale(:,:,1,i));
end



% save results to .gif file
for i_meas = 1:numMeas
for n = 1:cr
    rgbimage        =   zeros([size(im_ffd_scale,1),size(im_ffd_scale,2),3]);
    rgbimage(:,:,1) =   255*im_ffd_scale(:,:,i_meas,n);
    rgbimage(:,:,2) =   255*im_ffd_scale(:,:,i_meas,n);
    rgbimage(:,:,3) =   255*im_ffd_scale(:,:,i_meas,n);
    rgbimage        =   uint8(rgbimage);
    
    [imind,cm] = rgb2ind(rgbimage,256);
    % Write to the GIF File
    if (n == 1)&&(i_meas == 1)
        imwrite(imind,cm,[filename,'.tif'],'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,[filename,'.tif'],'gif','WriteMode','append');
    end
end
end

% save results to a video
numVideoFrame = cr*numMeas;
event_time = numMeas*cr*time_interval;   % time of the recorded event (s)
timeStamp = linspace(0,event_time,numVideoFrame);

% create video
hvid = VideoWriter([filename,'.avi']);
hvid.FrameRate = frame_rate;
open(hvid);

for i_meas = 1:numMeas
for n = 1:cr
    singleframe = im_ffd_scale(:,:,i_meas,n);
    singleframe_stamp = insertText(singleframe,[10 10],sprintf('%.2fms',timeStamp((i_meas-1)*cr+n)),'FontSize',30,'BoxColor',...
    [0 0 0],'BoxOpacity',0.4,'TextColor','white');
    writeVideo(hvid,singleframe_stamp);
end
end
close(hvid)



% save results to .mat file
recon = im_ffd_scale;
save(['recon_',filename,'.mat'],'recon');


end
