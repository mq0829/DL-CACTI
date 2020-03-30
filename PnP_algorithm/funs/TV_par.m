function savedata = TV_par(y, A,At,row,col,CSr,filename,P,Q,I)
close all
para.row = row;
para.col = col;
para.ori_im = I;
para.iter = 200;
para.lambda = 2;  
para.TVweight = 0.07; % the only parameter is this one, for challenge case, like the CSr<0.1, use a small one, as 0.1-0.2, for othter case, 0.1-1 is good

figure;

nr =1;
% now we use CSI_TV
im_IST    =   TV_CSI( y, nr, para, A, At);
%disp(['Running time: ' num2str(time)]);
PSNR_IST = psnr(im_IST,I);
subplot(2,3,1)
imshow(im_IST); axis off; title(['IST-TV, CSr: ' num2str(CSr)  ' PSNR: ' num2str(PSNR_IST)]);


nr =1;
% now we use GAP
para.iter = 200;
para.lambda = 1;  
im_GAP    =   TV_GAP( y, nr, para, A, At);

PSNR_GAP = psnr(im_GAP,I);
subplot(2,3,2)
imshow(im_GAP); axis off; title(['GAP-TV, CSr: ' num2str(CSr)  ' PSNR: ' num2str(PSNR_GAP)]);

% gap wavelet
ydim = round(CSr*row*col);


%%
stopc.iternum = 200;
stopc.err = 10^-5;
acc = 1;
block.row = 2;
block.col = 2;
m_star = ceil(ydim/(block.row*block.col));
[theta_wL21_dwt_tree model_wL21_dwt_tree] = GAP_2D_wL21_tree_newhadmard(y,  row,col, block,'wavelet', m_star,stopc,acc,P,Q);

theta_2d_wL21_dwt_tree = reshape(theta_wL21_dwt_tree,[row, col]);
im_GAP_dwt = theta_2d_wL21_dwt_tree/max(theta_2d_wL21_dwt_tree(:));
PSNR_gap_dwt = psnr(im_GAP_dwt,I);
subplot(2,3,3)
imshow(im_GAP_dwt); axis off; title(['GAP-wavelet, CSr: ' num2str(CSr)  ' PSNR: ' num2str(PSNR_gap_dwt)]);

%% TwiST
tau = 0.15;
% denoising function;
tv_iters = 5;
Psi = @(x,th)  tvdenoise(x,2/th,tv_iters);
% TV regularizer;
Phi = @(x) TVnorm(x);
tolA = 1e-4;
lam1=1e-4;  
AT = @(x) reshape(At(x), [row, col]);
A_twist = @(x) A(x(:));
% -- TwIST ---------------------------
% stop criterium:  the relative change in the objective function 
% falls below 'ToleranceA'
[im_twist,dummy,obj_twist,...
    times_twist,dummy,mse_twist]= ...
         TwIST(y,A_twist,tau,...
         'AT', AT, ...
         'lambda',lam1,...
         'True_x', I,...       
         'Psi', Psi, ...
         'Phi',Phi, ...
         'Monotone',1,...
         'StopCriterion',1,...
       	 'ToleranceA',tolA,...
         'Verbose', 0);
PSNR_twist = psnr(im_twist,I);
subplot(2,3,4)
imshow(im_twist); axis off; title(['TwIST-TV, CSr: ' num2str(CSr)  ' PSNR: ' num2str(PSNR_twist)]);

%
im_TVAL3 = run_TVAL3_nosave(A, At,y,row,col);
PSNR_TVAL3 = psnr(im_TVAL3,I);
subplot(2,3,5)
imshow(im_TVAL3); axis off; title(['TVAL3, CSr: ' num2str(CSr)  ' PSNR: ' num2str(PSNR_TVAL3)]);

savename = [filename '_TV_R' num2str(row) '_CSr' num2str(CSr)];
saveas(gcf,['.\Result\' savename '.fig']);
saveas(gcf,['.\Result\' savename '.png']);

savedata.im_IST = im_IST;
savedata.PSNR_IST = PSNR_IST;
savedata.im_GAP = im_GAP;
savedata.PSNR_GAP = PSNR_GAP;

savedata.im_GAP_dwt = im_GAP_dwt;
savedata.PSNR_gap_dwt = PSNR_gap_dwt;

savedata.im_twist = im_twist;
savedata.PSNR_twist = PSNR_twist;


savedata.im_TVAL3 = im_TVAL3;
savedata.PSNR_TVAL3 = PSNR_TVAL3;
% savedata.Im_TV_GAP = f_GAP_TV;
% savedata.psnr_TV_GAP = psnr_GAP_TV;




end