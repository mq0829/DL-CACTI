 function w = myfastidct2_givenT(theta,P2I,Bdct2I)
 
I_pat_3d_dct_shrink = Bdct2I(theta);
 im_out = P2I(I_pat_3d_dct_shrink);
 w  = im_out(:);
% w = reshape( shiftdim(T_col*shiftdim(T_row*reshape(theta_temp,[row col]),1),1),  [row*col 1]);

 end
