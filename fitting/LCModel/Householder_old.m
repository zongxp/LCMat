function [R,Q]=Householder_old(A)
% [R,Q]=Householder(A)
% A is m*n where m>n.  Q is the householder transformation matrix for upper triangular matrix.

Q=eye(size(A,1));
for j=1:size(A,2)
    [C,y,Qn]=H1(j,j+1,A(:,j),A(:,j+1:end));
   
    A(:,j)=y;
    A(:,j+1:end)=C;
    
    Q=Qn*Q;
end

R=A(1:size(A,2),:);



function [C,y,Q]=H1(p,l,nu,C)

m=length(nu);
s=nu(p)^2+sum(nu(l:m).^2);
s=sqrt(s);
if nu(p)>0
    s=-s;
end
if p>=l || p>m || p<1
    error('p, l, or m value error');
end


h=nu(p)-s;
nu(p)=s;
b=nu(p)*h;
y=zeros(size(nu));
y(1:l-1)=nu(1:l-1);

u=zeros(m,1);
u(p)=h;
u(l:m)=nu(l:m);

if b==0
    Q=eye(m);
else
    Q=eye(m)+u*u'/b;
end

if b==0 || ~exist('C','var') || isempty(C)
    return;
end

for j=1:size(C,2)
    s=(C(p,j)*h+sum(C(l:m,j).*nu(l:m)))/b;
    C(p,j)=C(p,j)+s*h;
    C(l:m,j)=C(l:m,j)+s*nu(l:m);
end



