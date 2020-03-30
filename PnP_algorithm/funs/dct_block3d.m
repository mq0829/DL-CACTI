function I_dct = dct_block3d(I_pat_3d,patch, T_dct_s, T_dct_t, cluster_indx,cluster,num_cluster)
%[patch, ~, patch_num] = size(I_pat_3d);

%nc = length(T_dct); % number of clusters
%if( cluster.eq )
%I_dct = zeros(size(I_pat_3d));
%else
I_dct = cell(1,num_cluster);
%end
for c=1:num_cluster
    I_pat_3d_dct_temp = I_pat_3d(:,:,cluster_indx{c});  % get the pathes in the same class 
    I_pat_3d_dct = reshape(T_dct_s*reshape(I_pat_3d_dct_temp, [patch, patch*cluster.cluster_num(c)]), [patch, patch, cluster.cluster_num(c)]);
    I_pat_3d_dct1 = reshape(T_dct_s*reshape(shiftdim(I_pat_3d_dct,1), [patch, cluster.cluster_num(c)*patch]), [patch,cluster.cluster_num(c), patch]);
    
    if(cluster.eq )
    I_pat_3d_dct2 = reshape(T_dct_t{c}*reshape(shiftdim(I_pat_3d_dct1,1), [cluster.cluster_num(c), patch*patch]), [cluster.cluster_num(c),patch, patch]);
    I_dct{c} = reshape(shiftdim(I_pat_3d_dct2,1), [patch*patch, cluster.cluster_num(c)]);
    else
    I_pat_3d_dct2 = reshape(shiftdim(I_pat_3d_dct1,1), [cluster.cluster_num(c), patch*patch]);
    I_pat_3d_dct3 = [I_pat_3d_dct2; zeros(cluster.cluster_num_ext(c) - cluster.cluster_num(c), patch^2)];
    I_dct{c} =  reshape(shiftdim(reshape(T_dct_t{c}*I_pat_3d_dct3, [cluster.cluster_num_ext(c), patch, patch]),1), [patch*patch, cluster.cluster_num_ext(c)]);
    end
end
end