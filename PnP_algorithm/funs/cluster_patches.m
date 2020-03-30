function [I2dct3d,dct2I3d] = cluster_patches(im, para,I2P,P2I,T_dct)
    row = para.row;
    col = para.col;
    patch = para.patch;
    step = para.step;
    
    f_3d =  image2patches3d(reshape(im,[row,col]), patch, patch, step, step);
          [cluster_indx,cluster_num,T_dct_t] = patch_cluster(f_3d,patch,para);
            if (strcmp(para.T_t,'dct'))
                cluster.eq = true;
                cluster.cluster_num = cluster_num;
                cluster.cluster_num_ext = cluster_num;
            else
                cluster.eq = false;
                cluster.cluster_num = cluster_num;
                cluster.cluster_num_ext = 2.^ceil(log2(cluster_num));
            end
           % now we get the idx for each patch
              I2Bdct3d = @(I) dct_block3d(I, patch, T_dct, T_dct_t, cluster_indx, cluster,para.cluster);
              Bdct2I3d = @(P) idct_block3d(P, patch, T_dct, T_dct_t, cluster_indx, cluster,para.cluster);
              I2dct3d = @(I) I2Bdct3d(I2P(I));
              dct2I3d = @(P) P2I(Bdct2I3d(P));
end