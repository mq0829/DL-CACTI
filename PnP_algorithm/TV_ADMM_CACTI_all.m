function  [im, PSNR_save]    =   TV_ADMM_CACTI_all( y, para, M_func,Mt_func)
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
    TVm = 'ITV3D_cham';
end

im          =    At( y );

eta = para.eta;
Phi_sum = para.Phi_sum;
%nc = size(im,3);
f      =   im;
b = 0;
theta=f;
%y1 = zeros(size(y));
for  iter = 1 : para.iter   
     if (mod(iter, 10) == 1)
        if isfield(para,'ori_im')
            PSNR     =   psnr( theta, para.ori_im);      
            PSNR_f     =   psnr( f, para.ori_im);      
            disp( [TVm ' ADMM Image Recovery, Iter ' num2str(iter) ': PSNR theta =' num2str(PSNR) ', PSNR f = ' num2str(PSNR_f)]); %, iter, PSNR, PSNR_f] );
        end
     end
   %  for ii = 1 : 1
         fpb = f+b;
            fb        =   A(fpb);
            theta = fpb + (At( (y-fb)./(Phi_sum+eta)));

     
        switch TVm
            case 'ATV_ClipA'
                f         =   TV_denoising(theta-b,  para.TVweight./eta,5);  
            case 'ATV_ClipB'
                f         =    TV_denoising_clip_LB(theta-b, para.TVweight./eta,5); 
                %f         =    TV_denoising_ClipB(theta-b,  para.TVweight,2e-2,50,5); 
            case 'ATV_cham'
                f         =     tvdenoise_cham_ATV2D(theta-b,  1/para.TVweight,5);  
            case 'ATV_FGP'
                f         =     fgp_denoise_ATV2D(theta-b,  para.TVweight,2);  
            case 'ITV2D_cham'
                f         =     tvdenoise_cham_ITV2D(theta-b,  1/para.TVweight,5);  
            case 'ITV2D_FGP'
                f         =     fgp_denoise_ITV2D(theta-b,  para.TVweight,2);  
            case 'ITV3D_cham'
                f         =     tvdenoise_cham_ITV3D(theta-b,  1/para.TVweight,5);  
            case 'ITV3D_FGP'
                f         =     fgp_denoise_ITV3D(theta-b,  para.TVweight,2);  
        end
        
        %end
        b = b -theta+f;
         PSNR_save(iter)     =   psnr( theta, para.ori_im); 
end
im     =  theta;

end