function x0 = TV_denoising_MM(y0,lambda,mu,iter)

if nargin<3
    iter = 20;
end

y_size = size(y0);

y_vec = y0(:);
%x_vec = zero(size(y_vec));

%z= zeros(max(size(y0)-1,1));  % the differential
%z= zeros(max(size(y0)-1,1));  % the differential
%alpha = 5;
m = y_size(1);
n = y_size(2);
T = y_size(3);
        % genera the TV matrix
Dv = TVmatrix(m,n,T,'V');
Dh = TVmatrix(m,n,T,'H');
Dt = TVmatrix(m,n,T,'T');
        
D = [Dv; Dh; mu*Dt];
DDT = sparse(D*D.');
Dy = D*y_vec;
z = zeros(3*m*n*T,1);

for it = 1:iter
    
    x_vec = y_vec - D.'*z;
     w = sqrt((Dv*x_vec).^2 + (Dh*x_vec).^2 + mu*(Dt*x_vec).^2);
    F = sparse(1:(3*m*n*T), 1:(3*m*n*T), [2*w./lambda; 2*w./lambda; 2*w./lambda]) + DDT;
    %x = y - D'*(F\Dy);
    %z = F\Dy;
    [z,~]=lsqr(@afun,Dy,1e-15,10,[],[],z); 

end
x0 = reshape(x_vec, [m,n,T]);
%end

function yv = afun(x,str)
    tempval= F*x;
     switch str
           case 'transp'
                 yv = tempval;
           case 'notransp'
                 yv = tempval;
      end
end
end

function opD=TVmatrix(m,n,T,str)

if str=='V' % This will give matrix for Horizontal Gradient
    D = spdiags([-ones(n,1) ones(n,1)],[0 1],n,n);
    D(n,:) = 0;
    D = kron(D,speye(m));
    D = kron(speye(T),D);
elseif str=='H' %This will give matrix for Verticle Gradient
   D = spdiags([-ones(m,1) ones(m,1)],[0 1],m,m);
   D(m,:) = 0;
   D = kron(speye(n),D);
   D = kron(speye(T),D);
elseif str=='T' %This will give matrix for Verticle Gradient
   D = spdiags([-ones(T,1) ones(T,1)],[0 1],T,T);
   D(T,:) = 0;
   D = kron(D,speye(m*n));
  % D = kron(speye(n),D);
end


opD=D;

end

function y = dv(x)
    y = diff(x);
end

function y = dh(x)
    y = diff(x,1,2);
end    
    
function y = dvt(x)
    y = [-x(1,:); -diff(x) ; x(end,:)];
end

function y = dht(x)
    y = [-x(:,1) -diff(x,1,2) x(:,end)];
end    
    
function y = dvt_3d(x)
    y = cat(1, -x(1,:,:), -diff(x) , x(end,:,:));
end

function y = dht_3d(x)
    y = cat(2, -x(:,1,:), -diff(x,1,2), x(:,end,:));
end    
function y = dvt_4d(x)
    y = cat(1, -x(1,:,:,:), -diff(x) , x(end,:,:,:));
end

function y = dht_4d(x)
    y = cat(2, -x(:,1,:,:), -diff(x,1,2), x(:,end,:,:));
        
end    
    
    
function y =clip(x, lambda)
    y = sign(x).*(min(abs(x),lambda));
end
