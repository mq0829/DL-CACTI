function  [im, PSNR_save]   =   TV_GAP_CACTI_FFDnetfast( y, para, M_func,Mt_func)
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

im          =    At( y );



Phi_sum = para.Phi_sum;
nc = size(im,3);
f      =   im;
TV_iter = para.TV_iter;

%b = 0;
%theta=f;
lambda      =   para.lambda;
y1 = zeros(size(y));
for isig = 1:(length(para.sigma)+1)
    if(isig==1)
   
    
    for iter = 1:TV_iter
     if (mod(iter, 10) == 1)
        if isfield(para,'ori_im')
            %PSNR     =   psnr( theta, para.ori_im);      
            PSNR_f     =   psnr( f, para.ori_im);      
            disp( [TVm ' GAP Image Recovery, Iter ' num2str(iter) ':, PSNR f = ' num2str(PSNR_f)]); %, iter, PSNR, PSNR_f] );
        end
     end
   %  for ii = 1 : 1
         fb        =   A( f );
            if(para.acc)
            y1 = y1+ (y-fb);
            f         =   f + lambda.*(At( (y1-fb)./Phi_sum ));
            else
                f         =   f + lambda.*(At( (y-fb)./Phi_sum ));
            end

     
        if(iter < TV_iter)
        switch TVm
            case 'ATV_ClipA'
                f         =   TV_denoising(f,  para.TVweight,5);  
            case 'ATV_ClipB'
                f         =    TV_denoising_clip_LB(f,  para.TVweight,5); 
               %  f         =    TV_denoising_ClipB(f,  para.TVweight,2e-2,50,5); 
            case 'ATV_cham'
                f         =     tvdenoise_cham_ATV2D(f,  1/para.TVweight,5);  
            case 'ATV_FGP'
                f         =     fgp_denoise_ATV2D(f,  para.TVweight,2);  
            case 'ITV2D_cham'
                f         =     tvdenoise_cham_ITV2D(f,  1/para.TVweight,5);  
            case 'ITV2D_FGP'
                f         =     fgp_denoise_ITV2D(f,  para.TVweight,2);  
            case 'ITV3D_cham'
                f         =     tvdenoise_cham_ITV3D(f,  1/para.TVweight,5);  
            case 'ITV3D_FGP'
                f         =     fgp_denoise_ITV3D(f,  para.TVweight,2);  
        end
        else
            
        end
        
        %end
        %b = b -theta+f;
         
    end
    else
     sig_use = para.sigma(isig-1);
     for iter = 1:para.maxiter(isig-1)
    % if (mod(iter, 10) == 1)
        if isfield(para,'ori_im')
            %PSNR     =   psnr( theta, para.ori_im);      
            PSNR_f     =   psnr( f, para.ori_im);      
            disp( [TVm ' GAP Image Recovery, Iter ' num2str(iter) ':, PSNR f = ' num2str(PSNR_f)]); %, iter, PSNR, PSNR_f] );
        end
    % end
   %  for ii = 1 : 1
         fb        =   A( f );
            if(para.acc)
            y1 = y1+ (y-fb);
            f         =   f + lambda.*(At( (y1-fb)./Phi_sum ));
            else
                f         =   f + lambda.*(At( (y-fb)./Phi_sum ));
            end
            
            for ncc = 1:nc
            f(:,:,ncc)         =  FFD_Net_Denoiserfast(squeeze(f(:,:,ncc)), sig_use,para.net); %(f,  para.TVweight,5);  TV_GAP_CACTI_FFDnet
            end
            %f = f./255;
            
     end
     PSNR_save(isig,iter)     =   psnr(f, para.ori_im);
    end
im     =  f;

end