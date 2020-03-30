function I_out = GMM_EVT(f,para)
   patch = para.patch;
   step = para.step;
   K = para.GMM_K;
   Rank = para.rank;
   row = para.row;
   col = para.col;


    I_pat = image2patches_fast(f, patch, patch, step, step);
    % now we train GMM
    options = statset('Display','off','MaxIter',200);
    obj = gmdistribution.fit(I_pat.',K,'Regularize',10^-3,'Options',options);
    pai = obj.PComponents; Mu = obj.mu'; Sig = obj.Sigma;

    % now we impose low rank of the Sig

    for kk = 1:K;
        sig_temp = Sig(:,:,kk);
        [V, D, W] = eig(sig_temp);
        d = diag(D);
        d = max(0, d-d(patch^2-Rank));
        D = diag(d);
%         V_low = V(:,end-Rank+1:end);
%         D_low = D(end-Rank+1:end,end-Rank+1:end);
        Sig_low(:,:,kk) = V*D*W';
    end

    % now we get the low rank of Sig, we need to compute the posterior of each
    % patch

    X = GMM_CS_Inv_samePhi(I_pat,eye(patch^2),Sig_low,Mu,pai);
    I_out = patches2image_fast(X, row,col, step, step);
end