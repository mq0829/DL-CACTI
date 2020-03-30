function [cluster_indx,cluster_num,T_dct_t] = patch_cluster(f_3d,patch,para)

           num_patch = size(f_3d,3);
           f_2d = reshape(f_3d,[patch^2, num_patch]);
           idx_patch = k_means(f_2d', 'random', para.cluster);
           for nc = 1:para.cluster
               cluster_indx{nc} = (idx_patch == nc);
               cluster_num(nc) = nnz(idx_patch == nc);
               if(cluster_num(nc)>0)
                   if (strcmp(para.T_t,'dct'))
                   T_dct_t{nc} = dctmtx(cluster_num(nc));
                   else
                       % this we use Haar wavelet
                      sig_level = ceil(log2(cluster_num(nc))); 
                      qmf   = MakeONFilter('Haar'); 
                      if(sig_level<8)
                          level = 3;
                      else
                          level = 5;
                      end
                      T_dct_t{nc} = get_waveletMatrix(qmf,sig_level,level,level);   
                   end
               end
           end
end