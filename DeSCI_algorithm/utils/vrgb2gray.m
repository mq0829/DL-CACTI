function vgray = vrgb2gray(vrgb)
%VRGB2GRAY Convert RGB video to gray.
nlast = size(vrgb,ndims(vrgb)); % size of the last dimension
odims = repmat({':'},1,ndims(vrgb)-1);
vgray = zeros(size(vrgb,1),size(vrgb,2),nlast);
for ilast = 1:nlast
    imrgb = vrgb(odims{:},ilast);
    vgray(:,:,ilast) = rgb2gray(imrgb);
end

end

