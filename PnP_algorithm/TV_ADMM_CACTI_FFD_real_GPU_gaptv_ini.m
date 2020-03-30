%%
function  im   =   TV_ADMM_CACTI_FFD_real_GPU( y, para, M_func,Mt_func)

% y = meas;
% M_func = A;
% Mt_func = At;

global recon_ini;

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
f      =   recon_ini;
b = zeros(size(f));
theta=f;
%y1 = zeros(size(y));
for isig = 1:(length(para.sigma))
    sig_use = para.sigma(isig);
    for iter = 1:para.maxiter(isig)
        if (mod(iter, 10) == 1)
            if isfield(para,'ori_im')
                PSNR     =   psnr( theta, para.ori_im);
                PSNR_f     =   psnr( f, para.ori_im);
                disp( [ 'FFDnet ADMM Image Recovery, Iter ' num2str(iter) ': PSNR theta =' num2str(PSNR) ', PSNR f = ' num2str(PSNR_f)]); %, iter, PSNR, PSNR_f] );
            else
                disp( [' FFDnet ADMM Image Recovery, Iter ' num2str(iter) ': L2 error =' num2str(norm(y-A(f))) ]);
            end
%             figure;imshow(theta(:,:,5)/max(theta(:)));
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
    
    % PSNR_save(iter)     =   psnr( theta, para.ori_im);
end
im     =  theta;

end