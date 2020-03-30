%get video data from a video file
%Input:
%   inname -- input video file name
%Ouput:
%   vData -- video data
function [vData]=read_video(inname)
video_d = VideoReader(inname);
nFrames = video_d.NumberOfFrames;
vidHeight = video_d.Height;
vidWidth = video_d.Width;
vData(1:nFrames) = struct('cdata', zeros(vidHeight, vidWidth, 3, 'uint8'), 'colormap', []);
 for k = 1 : nFrames
    vData(k).cdata = read(video_d, k);
 end
