function  [im, PSNR_save]    =   TV_TwIST_CACTI_all( y, para, M_func,Mt_func)
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

x         =    At( y );
%B  =size(x,3);


lam1=min(para.Phi_sum(:));   lamN=max(para.Phi_sum(:));

rho0 = (1-lam1/lamN)/(1+lam1/lamN);
alpha = 2/(1+sqrt(1-rho0^2));
beta  = alpha*2/(lam1+lamN);

xm2=x;
xm1=x;

for  iter = 1 : para.iter   
    if (mod(iter, 10) == 1)
        if isfield(para,'ori_im')
            PSNR     =   psnr(x, para.ori_im);                
            disp( [TVm ' TwIST Image Recovery, Iter ' num2str(iter) ':, PSNR Z = ' num2str(PSNR)]); %, iter, PSNR, PSNR_f] );
        end
    end
      z=x-At(A(x) - y);

     
        switch TVm
            case 'ATV_ClipA'
                 x         =   TV_denoising(z,  para.TVweight,5);  
            case 'ATV_ClipB'
                 x        =    TV_denoising_clip_LB(z,  para.TVweight,5); 
                 %x        =    TV_denoising_ClipB(z,  para.TVweight,2e-2,30,4); 
            case 'ATV_cham'
                 x         =     tvdenoise_cham_ATV2D(z,  1/para.TVweight,5);  
            case 'ATV_FGP'
                 x         =     fgp_denoise_ATV2D(z,  para.TVweight,2);  
            case 'ITV2D_cham'
                 x         =     tvdenoise_cham_ITV2D(z,  1/para.TVweight,5);  
            case 'ITV2D_FGP'
                 x         =     fgp_denoise_ITV2D(z,  para.TVweight,2);  
            case 'ITV3D_cham'
                 x         =     tvdenoise_cham_ITV3D(z,  1/para.TVweight,5);  
            case 'ITV3D_FGP'
                 x         =     fgp_denoise_ITV3D(z,  para.TVweight,2);  
        end
       x = (alpha-beta)*xm1 + (1-alpha)*xm2 + beta*x;
    %Z=X_iter+t_old/t_new*(Z_iter-X_iter)+(t_old-1)/t_new*(X_iter-X_old);
    xm2 = xm1;
    xm1 = x;
    PSNR_save(iter)     =   psnr(x, para.ori_im);  
end

im     =  x;

end