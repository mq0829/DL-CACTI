function x0 = TV_denoising_ClipB_debug(y0,lambda,rho,iter_out, iter_in, truth)

%rho = 1;% This is clipping vertically and horizontally and in the ADMM framework
v = 0;
u = 0;

if nargin<3
    iter = 100;
end
z= zeros(max(size(y0)-1,1));  % the differential
%z= zeros(max(size(y0)-1,1));  % the differential
alpha = 5;
for it = 1:iter_out
    if(size(y0,2)==1)  % 1D vector
    x0 = y0 - dvt(z);
    z = clip(z + 1/alpha*dv(x0), lambda/2);
    elseif(size(y0,2)>1 && (size(y0,3)==1))  % 2D image
        if(it ==1)
            zh = zeros(size(y0,1), size(y0,2)-1);
            zv = zeros(size(y0,1)-1, size(y0,2));
        end
    x0h = y0 - dht(zh);
    x0v = y0 - dvt(zv);
    x0 = (x0h + x0v)./2;
    zh = clip(zh + 1/alpha*dh(x0), lambda/2);
    zv = clip(zv + 1/alpha*dv(x0), lambda/2);
    elseif(size(y0,3)>1 && (size(y0,4)==1))  % 3D video or hyperspectral images, but ATV is done by each 2D frame. We foucs on thsi now
        if(it ==1)
            zh = zeros(size(y0,1), size(y0,2)-1, size(y0,3));
            zv = zeros(size(y0,1)-1, size(y0,2), size(y0,3));
        end
        
        for iit = 1:iter_in
        x0 = y0 + rho*(v-u)- dht_3d(zh);
        zh = clip(zh + 1/alpha*dh(x0), lambda/2);
        end
        for iit = 1:iter_in
        v = x0 + u - dvt_3d(zv);
        zv = clip(zv + 1/alpha*dv(v), lambda/rho/2);
        end
    
    u = u + (x0 -v);
    
   % rho = 0.9*rho;
    
    elseif(size(y0,4)>1 && (size(y0,5)==1))  % 4D hyperspectral-video, but TV is done by each 2D frame
        if(it ==1)
            zh = zeros(size(y0,1), size(y0,2)-1, size(y0,3) , size(y0,4));
            zv = zeros(size(y0,1)-1, size(y0,2), size(y0,3) , size(y0,4));
        end
    x0h = y0 - dht_4d(zh);
    x0v = y0 - dvt_4d(zv);
    x0 = (x0h + x0v)./2;
    zh = clip(zh + 1/alpha*dh(x0), lambda/2);
    zv = clip(zv + 1/alpha*dv(x0), lambda/2);
    end
   % x0 = x0./max(x0(:)+eps);
   
  disp(['Iter: ' num2str(it) ', PSNR: ' num2str(psnr(x0, truth)) ]);
end





function y = dv(x)
    y = diff(x);


function y = dh(x)
    y = diff(x,1,2);
    
    
function y = dvt(x)
    y = [-x(1,:); -diff(x) ; x(end,:)];


function y = dht(x)
    y = [-x(:,1) -diff(x,1,2) x(:,end)];
    
    
function y = dvt_3d(x)
    y = cat(1, -x(1,:,:), -diff(x) , x(end,:,:));


function y = dht_3d(x)
    y = cat(2, -x(:,1,:), -diff(x,1,2), x(:,end,:));
    
function y = dvt_4d(x)
    y = cat(1, -x(1,:,:,:), -diff(x) , x(end,:,:,:));


function y = dht_4d(x)
    y = cat(2, -x(:,1,:,:), -diff(x,1,2), x(:,end,:,:));
        
    
    
    
function y =clip(x, lambda)
    y = sign(x).*(min(abs(x),lambda));
