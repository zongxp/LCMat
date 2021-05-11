function [Uc,Ua,b_H0]=get_commute(U,H)
%[Uc,Ua,b_H0]=get_commute(U,H)
% calculate the component of U that commute or anticommute with H
% b_H0: column eigenvectors for zero eigenvalues of H.
[v,d]=eig(H);

U2=v'*U*v;
H2=v'*H*v;

ia=[];
ib=[];

for i=1:size(d,1)
    if d(i,i)~=0
        ia(end+1)=i;
    else
        ib(end+1)=i;
    end
end

n0=length(ia);
Usub=zeros(n0,n0);
Hsub=zeros(n0,n0);

for i=1:length(ia)
     Usub(i,:)=U2(ia(i),ia);
     Hsub(i,:)=H2(ia(i),ia);
end

Usub_c=Usub/2+1/2*Hsub*Usub/Hsub;

Uc=zeros(size(d));
Uc(ia,ia)=Usub_c;
Uc(ib,ib)=U2(ib,ib);

Ua=U-Uc;

Ua=v*Ua*v';
Uc=v*Uc*v';

b_H0=v(:,ib);


    
    

