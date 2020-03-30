function [T,T_all] = get_waveletMatrix(h,sig_level,level,decompose_level)
% T = get_waveletMatrix(sig_level,level,decompose_level)
% construct the matrix of wavelet tranform
% Copyright (C) 2008-2012 by Xuejun Liao @ Duke University
% Comments to:  xuejun.liao@gmail.com, xjliao@ee.duke.edu, (919)-660-5548
% 11 November 2008

%hgcode=4; h=choosehg(hgcode);

%sig_level = 4;
%level = 3;

g = zeros(size(h)); 
L = length(h);
for i=1:L
    g(i)=(-1)^(i-1)*h(L-i+1); 
end

T = eye(2^sig_level);
T_all = {T};
for t=sig_level:-1:max(level,decompose_level+1)
  T = formhg(0,sig_level,t,h,g)*T;
  T_all = [T_all {formhg(0,sig_level,t,h,g)}];
end

