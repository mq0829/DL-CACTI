function  [im, PSNR_save]    =   TV_FISTA_CACTI_all( y, para, M_func,Mt_func)
if nargin==4%function
    A=@(x) M_func(x);
    At=@(z) Mt_func(z);
else%Matrix
    A=@(x)M_func*x;
    At=@(z)M_func'*z;
end


if (isfield(para,'TVm'))
    TVm = para.TVm;
else
    TVm = 'ATV_ClipA';  % this is the fastest algorithm and usually gives a good result
end

X_iter         =    At( y );
B  =size(X_iter,3);
L = B;

Z = X_iter;
t_new=1; 
for  iter = 1 : para.iter   
    if (mod(iter, 10) == 1)
        if isfield(para,'ori_im')
            PSNR     =   psnr(Z, para.ori_im);                
            disp( [TVm ' FISTA Image Recovery, Iter ' num2str(iter) ':, PSNR Z = ' num2str(PSNR)]); %, iter, PSNR, PSNR_f] );
        end
     end
    X_old=X_iter;
    t_old=t_new;
    % gradient step
   
    Z=Z-1/L*At(A(Z) - y);

     
        switch TVm
            case 'ATV_ClipA'
                 X_iter         =   TV_denoising(Z,  para.TVweight,5);  
            case 'ATV_ClipB'
                 X_iter         =    TV_denoising_clip_LB(Z,  para.TVweight,5); 
                 %X_iter         =    TV_denoising_ClipB(Z,  para.TVweight,2e-2,30,4); 
            case 'ATV_cham'
                 X_iter         =     tvdenoise_cham_ATV2D(Z,  1/para.TVweight,5);  
            case 'ATV_FGP'
                 X_iter         =     fgp_denoise_ATV2D(Z,  para.TVweight,2);  
            case 'ITV2D_cham'
                 X_iter         =     tvdenoise_cham_ITV2D(Z,  1/para.TVweight,5);  
            case 'ITV2D_FGP'
                 X_iter         =     fgp_denoise_ITV2D(Z,  para.TVweight,2);  
            case 'ITV3D_cham'
                 X_iter         =     tvdenoise_cham_ITV3D(Z,  1/para.TVweight,5);  
            case 'ITV3D_FGP'
                 X_iter         =     fgp_denoise_ITV3D(Z,  para.TVweight,2);  
        end
      t_new=(1+sqrt(1+4*t_old^2))/2;
    %Z=X_iter+t_old/t_new*(Z_iter-X_iter)+(t_old-1)/t_new*(X_iter-X_old);
     Z=X_iter+(t_old-1)/t_new*(X_iter-X_old);
     PSNR_save(iter) = psnr(Z, para.ori_im);
end

im     =  Z;

end