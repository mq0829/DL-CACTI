function I_dct = dct_block2d(I_pat_3d,T_dct)
[patch, ~, patch_num] = size(I_pat_3d);
I_pat_3d_dct = reshape(T_dct*reshape(I_pat_3d, [patch, patch*patch_num]), [patch, patch, patch_num]);
   I_pat_3d_dct1 = reshape(T_dct*reshape(shiftdim(I_pat_3d_dct,1), [patch, patch_num*patch]), [patch, patch_num, patch]);
   I_dct = shiftdim(I_pat_3d_dct1,2);
end