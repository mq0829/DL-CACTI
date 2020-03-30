function [PSNR, PSNR_mos] = my_genneratemoive(img, Row, Col, CodeFrame, ColT,X,  y,X_ori_color,X_mos,fname,Algorim)
Img_recon_sensor = zeros(Row,Col,CodeFrame*ColT);
   
for nR = 1:(Row/2)
    for nC = 1:(Col/2)
        Img_recon_sensor((nR-1)*2+1,(nC-1)*2+1,:) = img(nR,nC,:,1);
        Img_recon_sensor((nR-1)*2+1,(nC-1)*2+2,:) = img(nR,nC,:,2);
        Img_recon_sensor((nR-1)*2+2,(nC-1)*2+1,:) = img(nR,nC,:,3);
        Img_recon_sensor((nR-1)*2+2,(nC-1)*2+2,:) = img(nR,nC,:,4);
    end
end

for k=1:(CodeFrame*ColT)
     temp = uint8(Img_recon_sensor(:,:,k)/max(max(Img_recon_sensor(:,:,k)))*255);
     X_recon_col(:,:,:,k) = demosaic(temp,'rggb');
     PSNR(k) = psnr(X_ori_color(:,:,:,k),im2double(X_recon_col(:,:,:,k)));
     PSNR_mos(k) = psnr(X(:,:,k),Img_recon_sensor(:,:,k));
 end

figure; 
plot( PSNR_mos,'r-*','LineWidth',4,'MarkerSize',10);
hold on;
plot(PSNR,'k->','LineWidth',4,'MarkerSize',10);

xlabel('Frame','fontsize',20);
ylabel('PSNR (dB)','fontsize',20);
set(gca,'fontsize',20);
legend(['Reconstructed mosaic video, average PSNR:' num2str(mean(PSNR_mos))],['Reconstructed demosaic video, average PSNR:' num2str(mean(PSNR))]);

%  figure;
%  for k=1:(CodeFrame*ColT)
%      imshow(cat(2,X_ori_color(:,:,:,k),im2double(X_recon_col(:,:,:,k)))); title(['Frame: ' num2str(k) '  PSNR:' num2str(PSNR(k))]); pause(0.04);
%  end

% 
 savename = [fname '_T' num2str(ColT) '_F' num2str(CodeFrame) '_' Algorim];
 format short
 save(savename, '-v7.3', 'X_recon_col','Img_recon_sensor','PSNR','PSNR_mos');
 
 writerObj = VideoWriter([savename '.mp4'],'MPEG-4');
%writerObj.FrameRate = round(T/7);
writerObj.FrameRate = 12;

open(writerObj);
%scrsz = get(0,'ScreenSize');
%fig=figure('Position',[50 100 floor(scrsz(3)*0.8) floor(scrsz(4)*0.6)]);
%fig=figure;
fig=figure('Position',[50 100 1920 1088]);
for nF=1:(CodeFrame*ColT)
    

     subplot(2,3,1);
    imshow(X_ori_color(:,:,:,nF)); 
    title(['Original Video, Frame:' num2str(nF)]);
     freezeColors;
    
    subplot(2,3,2);
    imshow(X_mos(:,:,nF)/max(max(X_mos(:,:,nF)))); 
    title('Original mosaic video');
     freezeColors;
    
    subplot(2,3,3);
    imshow(y(:,:,ceil(nF/ColT))/max(max(y(:,:,ceil(nF/ColT))))); 
    title(['Measurement mosaic video, Frame: ' num2str(ceil(nF/ColT))]);
     freezeColors;
    
    subplot(2,3,4);
    imshow(X_recon_col(:,:,:,nF)); 
    title(['Recon demosaic video, PSNR: ' num2str(PSNR(nF))]);
     freezeColors;
     
     subplot(2,3,5);
    imshow(Img_recon_sensor(:,:,nF)/max(max(Img_recon_sensor(:,:,nF)))); 
    title(['Recon mosaic video, PSNR: ' num2str(PSNR_mos(nF))]);
     freezeColors;
     
    %pause(0.02); 
    frame = getframe(fig);
    writeVideo(writerObj,frame);
end
close(writerObj);
