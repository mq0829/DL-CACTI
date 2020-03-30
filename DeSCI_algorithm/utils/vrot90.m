function vrot = vrot90( v, rotnum )
%VROT90 Rotate the video of multiple 90 degrees.
if nargin < 2
    rotnum = 1;
end
if ndims(v) <= 3
    nframe = size(v,3);
    for iframe = 1:nframe
        vrot(:,:,iframe) = rot90(v(:,:,iframe),rotnum);
    end
else
    nframe = size(v,4);
    ncolor = size(v,3);
    for iframe = 1:nframe
        for icolor = 1:ncolor
            vrot(:,:,icolor,frame) = rot90(v(:,:,icolor,iframe),rotnum);
        end
    end
end
    
    

