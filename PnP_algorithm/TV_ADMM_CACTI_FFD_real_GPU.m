%%
function  im   =   TV_ADMM_CACTI_FFD_real_GPU( y, para, M_func,Mt_func)
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
nc = size(im,3);

TV_iter = para.TV_iter;

eta = para.eta;
Phi_sum = para.Phi_sum;
%nc = size(im,3);
f      =   im;
b = 0;
theta=f;
%y1 = zeros(size(y));
for isig = 1:(length(para.sigma)+1)
    if(isig==1)
     for iter = 1:TV_iter
     if (mod(iter, 10) == 1)
        if isfield(para,'ori_im')
            PSNR     =   psnr( theta, para.ori_im);      
            PSNR_f     =   psnr( f, para.ori_im);      
            disp( [TVm ' ADMM Image Recovery, Iter ' num2str(iter) ': PSNR theta =' num2str(PSNR) ', PSNR f = ' num2str(PSNR_f)]); %, iter, PSNR, PSNR_f] );
        else
            disp( [TVm ' ADMM Image Recovery, Iter ' num2str(iter) ': L2 error =' num2str(norm(y-A(f))) ]);
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
        
     end
    
     
     
     else
     sig_use = para.sigma(isig-1);
     for iter = 1:para.maxiter(isig-1)
        if (mod(iter, 10) == 1)
            if isfield(para,'ori_im')
                PSNR     =   psnr( theta, para.ori_im);      
                PSNR_f     =   psnr( f, para.ori_im);      
                disp( [ 'FFDnet ADMM Image Recovery, Iter ' num2str(iter) ': PSNR theta =' num2str(PSNR) ', PSNR f = ' num2str(PSNR_f)]); %, iter, PSNR, PSNR_f] );
            else
                disp( [' FFDnet ADMM Image Recovery, Iter ' num2str(iter) ': L2 error =' num2str(norm(y-A(f))) ]);
            end
            
         end
       %  for ii = 1 : 1
               fpb = f+b;
                fb        =   A(fpb);
                theta = fpb + (At( (y-fb)./(Phi_sum+eta)));

                for ncc = 1:nc
                f(:,:,ncc)         =  FFD_Net_DenoiserGPU(squeeze(theta(:,:,ncc)-b(:,:,ncc)), sig_use,para.net); %(f,  para.TVweight,5);  TV_GAP_CACTI_FFDnet
                end
                %f = f./255;
               b = b -theta+f;  
      end
   %  PSNR_save(isig,iter)     =   psnr(f, para.ori_im);
    end
        % PSNR_save(iter)     =   psnr( theta, para.ori_im); 
end
im     =  theta;

end