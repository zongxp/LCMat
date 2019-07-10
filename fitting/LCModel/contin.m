function [x,cov,dof,V]=contin(A,y,R,r,D,d,alpha)
% [x,cov,dof,V]=contin(A,y,R,r,D,d,alpha)
% A:  the design matrix, Ny*Nx
% y: the data, Ny*1
% R: the regularizor matrix, Nreg*Nx
% r: regularizor vector, Nreg*1
% D: inequality matrix, Nieq*Nx; Dx>=d
% d: inequality vector, Nieq*1
% alpha: regularizor parameter
% x: solution
% sigma: covariance matrix with unit noise
% dof: effective degree of freedom.
% V: residual norm
[C,HA]=Householder(A);

Nx=size(A,2);
Ny=size(A,1);
if Nx>Ny || size(A,2) ~= size(R,2)
    error('dimension error');
end

eta=HA(1:Nx,:)*y;

Nxe=Nx;
Nreg=size(R,1);
if Nreg<Nxe
    R(end+1:Nxe,:)=0;
    r(end+1:Nxe)=0;
elseif Nreg>Nxe
    error('Do not support Nreg>Nex yet');
end

[K2,H1,Z,W,S,U,Q] = remove_eqcon(C,R);

s=diag(S);
St=diag(sqrt(s.^2+alpha^2));
r1=U'*r;
gamma=Q'*(eta-C*Z*inv(H1)*r1);
gammat=gamma.*s./sqrt(s.^2+alpha^2); %gamma_tilde
option=optimset('LargeScale','off','Display','off');
xi=lsqlin(eye(Nxe),zeros(Nxe,1),-D*Z*inv(H1)*W*inv(St),-d+D*Z*inv(H1)*(r1+W*inv(St)*gammat),[],[],[],[],[],option);
%xi=lsqlin(eye(Nxe),zeros(Nxe,1),-D*Z*inv(H1)*W*inv(St),-d+D*Z*inv(H1)*(r1+W*inv(St)*gammat));

x=Z*inv(H1)*(W*inv(St)*(xi+gammat)+r1);

%% find error in x
binding=find(abs(D*x-d)<1e-6);

if ~isempty(binding)
    E = D(binding,:);
    [K2,H1,Z,W,S] = remove_eqcon(C,R,E);
end

s=diag(S);
G=diag(s./(s.^2+alpha^2));
dof=sum(s.^2./(s.^2+alpha^2));
V = sum((y-A*x).^2)+alpha^2*sum((r-R*x).^2);

tmp=(K2*Z*inv(H1)*W);
cov=tmp*G*G*tmp';

function [K2,H1,Z,W,S,U,Q] = remove_eqcon(C,R,E)

Nx=size(C,2);
if exist('E','var')
 Neq=size(E,1);
 [tmp,K]=Householder(E');
 K=K';
 K2=K(:,Neq+1:end);
 Nxe=Nx-Neq;
else
 Nxe=Nx;    
 K2=eye(Nxe);  %no equality constrains.
end

[U,H,Z]=svd(R*K2);
H1=H(1:Nxe,1:Nxe);
H1d=diag(H1);
H1d(H1d<1e-6*max(H1d))=1e-6*max(H1d);
H1=diag(H1d);
[Q,S,W]=svd(C*K2*Z*inv(H1));
