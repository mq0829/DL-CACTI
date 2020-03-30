function  [im, data_save]    =   TV_GAP_CACTI_cham3d( y, nr,para, M_func,Mt_func)
if nargin==5%function
    A=@(x) M_func(x);
    At=@(z) Mt_func(z);
else%Matrix
    A=@(x)M_func*x;
    At=@(z)M_func'*z;
end
%row = para.row;
%col = para.col;

im          =    At( y );
mu = para.mu;
%im          =    reshape(im,[row col]);
lamada      =   para.lambda;
%nc = size(im,3);
f      =   im;
y1 = zeros(size(y));
for  iter = 1 : para.iter   
     if (mod(iter, 10) == 1)
        if isfield(para,'ori_im')
            % PSNR     =   psnr( f./max(f(:)), para.ori_im./max(para.ori_im(:))); 
            PSNR     =   psnr( f, para.ori_im); 
            fprintf( 'TVcham3D-GAP Compressive Image Recovery, Iter %d : PSNR = %f\n', iter, PSNR );
        end
     end
     for ii = 1 : 1
            fb        =   A( f );
            if(para.acc)
            y1 = y1+ (y-fb);
            f         =   f + lamada.*(At( (y1-fb)./para.Phi_sum ));
            else
                f         =   f + lamada.*(At( (y-fb)./para.Phi_sum ));
            end
     end         
     
        %for nn = 1:nc 
        % f         =   tvdenoise_cham3d(f,  1./para.TVweight,mu,4);  
       %   f         =   tvdenoise_chammax_L13D(f,  1./para.TVweight,mu,4); 
        %  f         =   tvdenoise_chammax_L13D(f,  1./para.TVweight,4); 
          f         =   tvdenoise_chammax_L1w3D(f,  1./para.TVweight,mu,4); 
        %end
        data_save.psnr(iter) = psnr( f, para.ori_im); 
end
im     =  f;

end