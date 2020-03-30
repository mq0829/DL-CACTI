function  im_out   =  DCT_threshold_3d_block( I, CSr,I2dct, dct2I,num_cluster)

    I_pat_3d_dct = I2dct(I);
   
   % now we do theoresholding on each cluster
   for nc = 1:num_cluster
       temp = I_pat_3d_dct{nc};
       [I_pat_dct_sort, ind] = sort(abs(temp(:)), 'descend');   
       %CSr = 0.2;
       Re_num = round(CSr*length(temp(:)));

       if(Re_num < length(temp(:)))
       thero = I_pat_dct_sort(Re_num+1);
       else
           thero = 0;
       end

       temp = sign(temp).*max(abs(temp)-thero,0);
       I_pat_dct_shrink_cell{nc} = temp;
   end

  im_out = dct2I(I_pat_dct_shrink_cell);

