function [X, Recon] = GMM_CS_Inv_samePhi(y,Phi,Sig,mu,weight)
% here y is a matrix
     [m M] = size(Phi);  % Size of the Projection Matrix
     N = size(y,2);
      K = size(mu,2); % GMM component number
      Rinv = 1e-4*eye(m); 
      % Compute the distribution  
      %py = zeros(1,K);
      % Compute the mean and convariance
      for k=1:K

       P1 = inv(Phi*Sig(:,:,k)*Phi'+Rinv);
       P = (P1+P1')/2;
       %Recon.mu(:,:,k) = (Sig(:,:,k)*Phi')*P*(y-Phi*mu(:,k))+mu(:,k);
       Recon.mu(:,:,k) = bsxfun(@plus,(Sig(:,:,k)*Phi')*P*bsxfun(@minus,y,Phi*mu(:,k)),mu(:,k));
       res = bsxfun(@minus, y ,Phi*mu(:,k));
       likeli = -0.5*m*log(2*pi) + sum(log(diag(chol(P)))) - 0.5*sum(res.*(P*res),1);   
       logpy(:,k) = likeli + log(weight(k)+eps);

      end
      s = logsumexp(logpy');
      pymin = bsxfun(@minus,logpy',s);
      % sumpy = sum(py);
     
      Recon.weight = reshape(exp(pymin(:)),[K,N]); 
      
      
%       for k=1:K
%           EstX(:,k) = Recon.weight(k)*Recon.mu(:,k);
%       end
      Recon.X = squeeze(sum(bsxfun(@times, Recon.weight', shiftdim(Recon.mu,1)),2));
      X = Recon.X';
      
end