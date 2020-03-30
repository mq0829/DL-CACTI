function T=formhg(chpack,level,t,h,g)
% FORMHG: T=formhg(chpack,level,t,h,g)
%         Form the transform matrix T from h and g, refer for further
%         details to wavede.m, by which this sub routine is called.
%         chpack: 1 chooses wavelet packet decomposion
%                 0 chooses normal wavelet decomposion
%         level:  the length in power of 2 of input signal 
%         t:      current decomposition level
%         h:      low-pass filter coefficients
%         g:      high-pass filter coefficients
%         T:      transform matrix 
% Copyright (c) by Xuejun Liao
% December 1997

H=zeros(2^t,2^t/2);  G=H;
for i=1:2^t/2  % form the low-pass H and high-pass filter G  
  L=length(h);
  if 2^t-2*(i-1)>=L
    H((1+2*(i-1)):(L+2*(i-1)),i)=h;
    G((1+2*(i-1)):(L+2*(i-1)),i)=g;
  else
    shift=L-(2^t-2*(i-1));
    H((1+2*(i-1)):(L-shift+2*(i-1)),i)=h(1:(L-shift));
    H(1:shift,i)=h((L-shift+1):L);
    G((1+2*(i-1)):(L-shift+2*(i-1)),i)=g(1:(L-shift));
    G(1:shift,i)=g((L-shift+1):L);
  end
end
T=zeros(2^level);
% contatenate H to G so as to form one transform matrix T.
% In the wavelet case (choice==1), the upper left sub principal
% matrix is [H';G'], and the complement sub principal-diagonal
% matrix is left to be an anidentity matrix. In the wavelet packet
% case (choice==2), [H';G'] is lined recursively along the principal-
% diagonal.
if chpack==0
  T(1:2^t/2,1:2^t)=H';
  T((2^t/2+1):2^t,1:2^t)=G';
  for i=(2^t+1):(2^level)  
    T(i,i)=1;	           
  end
else
  for i=0:(2^(level-t)-1)
    T((2^t*i+1):(2^t*i+2^t/2),(2^t*i+1):(2^t*i+2^t))=H';
    T((2^t*i+2^t/2+1):(2^t*i+2^t),(2^t*i+1):(2^t*i+2^t))=G';
  end
end

