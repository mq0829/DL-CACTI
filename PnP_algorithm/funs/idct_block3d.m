function I_dct = idct_block3d(I_pat_3d,patch, T_dct_s, T_dct_t, cluster_indx,cluster,num_cluster)
%[patch, ~, patch_num] = size(I_pat_3d);

%nc = length(T_dct); % number of clusters

I_dct = zeros(patch, patch, sum(cluster.cluster_num));

for c=1:num_cluster
    I_pat_3d_dct_temp = reshape(I_pat_3d{c}, [patch, patch, cluster.cluster_num_ext(c)]);  % get the pathes in the same class 
    I_pat_3d_dct = reshape(T_dct_s'*reshape(I_pat_3d_dct_temp, [patch, patch*cluster.cluster_num_ext(c)]), [patch, patch, cluster.cluster_num_ext(c)]);
    I_pat_3d_dct1 = reshape(T_dct_s'*reshape(shiftdim(I_pat_3d_dct,1), [patch, cluster.cluster_num_ext(c)*patch]), [patch,cluster.cluster_num_ext(c), patch]);
    
    if(cluster.eq )
    I_pat_3d_dct2 = reshape(T_dct_t{c}'*reshape(shiftdim(I_pat_3d_dct1,1), [cluster.cluster_num_ext(c), patch*patch]), [cluster.cluster_num_ext(c),patch, patch]);
    I_dct(:,:,cluster_indx{c}) = shiftdim(I_pat_3d_dct2,1);
    else
    I_pat_3d_dct2 = reshape(shiftdim(I_pat_3d_dct1,1), [cluster.cluster_num_ext(c), patch*patch]);
    temp =  shiftdim(reshape(T_dct_t{c}'*I_pat_3d_dct2, [cluster.cluster_num_ext(c), patch, patch]),1);
    temp1 =  temp(:,:,1:cluster.cluster_num(c));
    I_dct(:,:,cluster_indx{c}) = temp1;
    end
end
end