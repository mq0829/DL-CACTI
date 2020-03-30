

function [im,  M_weight]= patches3d2image(M_patches, N1, N2, delta1, delta2)
% Reput patches in images in the appropriate positions with appropriate
% weight. It is the inverse procedure of image2patches3d
% M_patches is 3D with size patch x patch x num if patches
% Xin Yuan, 2015-05-01

if ( N1 > 65535 ) || ( N2 > 65535 ) 
    error('The image size should be smaller than 65535 x 65535');
end

[n1, n2, num_patches] = size(M_patches);



% % % n1 = 110;
% % % n2 = 1;

% the coordinates of the top-left point in all the patches are computed and
% stored in (XstartI, YstartI). XstartI or YstartI is a vector of length
% #patches. 
Xstart = uint16(1 : delta1 : N1 - n1 + 1);
Xstart = [Xstart  Xstart(end)+1:(N1 - n1 + 1)];
Ystart = uint16(1 : delta2 : N2 - n2 + 1);
Ystart = [Ystart  Ystart(end)+1:(N2 - n2 + 1)];
[YstartI, XstartI] = meshgrid(Ystart, Xstart);
YstartI = YstartI(:);
XstartI = XstartI(:);

n1_minus1 = n1 - 1;
n2_minus1 = n2 - 1;


im = zeros(N1, N2);
M_weight = zeros(N1, N2);

% use (one-layer) loop to extract the patches. This loop is inevitable in
% the reconstruction phase (patches2image) because we need to add the
% patches and accumulate the weight. 
for k = 1 : num_patches
    coor_x = XstartI(k):XstartI(k)+n1_minus1;
    coor_y = YstartI(k):YstartI(k)+n2_minus1;
    im(coor_x, coor_y) = im(coor_x, coor_y) + M_patches(:, :,k);
    M_weight(coor_x, coor_y) = M_weight(coor_x, coor_y) + 1;
end

im = im ./ M_weight;

