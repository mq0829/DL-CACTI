function [rim] = imcropmag(im,opts)
%IMCROPMAG Crop image with magnified image along the side.
if nargin<2
    opts = [];
end
[h,w,~] = size(im);
% [0] set default value of cropped image and the magnified image
crect   = round([h/2-h/20,w/2-w/20,h/10,w/10]); % rectangle position of the cropped image [h1 w1 height width]
magpos  = 'downright'; % position to put the magnified image "downright" for default
magsize = round([h/3 w/3]); % size of the magnified image [height width]
clinewidth = 2; % line width of the cropped image
mlinewidth = 4; % line width of the magnified image
% assign manual options
if isfield(opts,'crect')           crect = opts.crect;      end
if isfield(opts,'magpos')         magpos = opts.magpos;     end
if isfield(opts,'magsize')       magsize = opts.magsize;    end
if isfield(opts,'clinewidth') clinewidth = opts.clinewidth; end
if isfield(opts,'mlinewidth') mlinewidth = opts.mlinewidth; end

% [1] get the cropped image and have the magnified one displayed on the
% original image
cim = im(crect(1):crect(1)+crect(3)-1,crect(2):crect(2)+crect(4)-1,:); % crop image
mim = imresize(cim,magsize,'nearest');

rim = insertShape(im,'rectangle',crect([2 1 4 3]),'LineWidth',clinewidth,'Color','red');

exim = zeros([magsize+2*mlinewidth size(im,3)],'like',im);
exim(:,:,1) = 255; exim(:,:,2) = 0; exim(:,:,3) = 0; % RGB channel
if size(mim,3)==1
    mim = repmat(mim,[1 1 3]);
end
exim(1+mlinewidth:magsize(1)+mlinewidth,1+mlinewidth:magsize(2)+mlinewidth,:) = mim;

switch magpos % position of the magnified image
    case 'upleft'
        rim(1:magsize(1)+2*mlinewidth,1:magsize(2)+2*mlinewidth,:) = exim;
    case 'upright'
        rim(1:magsize(1)+2*mlinewidth,end-magsize(2)-2*mlinewidth+1:end,:) = exim;
    case 'downleft'
        rim(end-magsize(1)-2*mlinewidth+1:end,1:magsize(2)+2*mlinewidth,:) = exim;
    case 'downright'
        rim(end-magsize(1)-2*mlinewidth+1:end,end-magsize(2)-2*mlinewidth+1:end,:) = exim;
    otherwise
        error('Unsupported position of the magnified image!');
end

end

