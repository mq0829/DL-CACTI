function u = tvdenoise_oct_sp(f,lambda,iters)
%TVDENOISE  Total variation grayscale and color image denoising
%   u = TVDENOISE(f,lambda) denoises the input image f.  The smaller
%   the parameter lambda, the stronger the denoising.
%
%   The output u approximately minimizes the Rudin-Osher-Fatemi (ROF)
%   denoising model
%
%       Min  TV(u) + lambda/2 || f - u ||^2_2,
%        u
%
%   where TV(u) is the total variation of u.  If f is a color image (or any
%   array where size(f,3) > 1), the vectorial TV model is used,
%
%       Min  VTV(u) + lambda/2 || f - u ||^2_2.
%        u
%
%   TVDENOISE(...,Tol) specifies the stopping tolerance (default 1e-2).
%
%   The minimization is solved using Chambolle's method,
%      A. Chambolle, "An Algorithm for Total Variation Minimization and
%      Applications," J. Math. Imaging and Vision 20 (1-2): 89-97, 2004.
%   When f is a color image, the minimization is solved by a generalization
%   of Chambolle's method,
%      X. Bresson and T.F. Chan,  "Fast Minimization of the Vectorial Total
%      Variation Norm and Applications to Color Image Processing", UCLA CAM
%      Report 07-25.
%
%   Example:
%   f = double(imread('barbara-color.png'))/255;
%   f = f + randn(size(f))*16/255;
%   u = tvdenoise(f,12);
%   subplot(1,2,1); imshow(f); title Input
%   subplot(1,2,2); imshow(u); title Denoised

% Pascal Getreuer 2007-2008
%  Modified by Jose Bioucas-Dias  & Mario Figueiredo 2010 
%  (stopping rule: iters)
%   


dt = 0.25;
f_s = single(f);
N = size(f);
id = [2:N(1),N(1)];
iu = [1,1:N(1)-1];
ir = [2:N(2),N(2)];
il = [1,1:N(2)-1];
p1 = zeros(size(f), 'single');
p2 = zeros(size(f), 'single');
divp = zeros(size(f), 'single');



%     
     for i=1:iters
        z = divp - f_s*single(lambda);
        z1 = z(:,ir,:) - z;
        %z1 = cat(2, diff(z,1,2), zeros(N(1), 1, N(3)));
        z2 = z(id,:,:) - z;
        %z2 = cat(1, diff(z,1,1), zeros(1, N(2), N(3)));
        denom = 1 + dt*sqrt(sum(z1.^2 + z2.^2,3));
        p1 = bsxfun(@rdivide, (p1 + dt*z1),denom);
        p2 = bsxfun(@rdivide, (p2 + dt*z2),denom);
        divp = p1 - p1(:,il,:) + p2 - p2(iu,:,:);
    end



u = f - divp/single(lambda);
