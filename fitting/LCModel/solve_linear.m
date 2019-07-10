function [V,c_b,ec_b,p,dof]=solve_linear(sp,model_base,bspl,Rb,alpha,sigma)
%[c_b,ec_b,p,V]=solve_linear(sp,model_base,bspl,Rb,alpha,sigma)

Nm=size(model_base,2);
Nb=size(bspl,2);
A=[model_base,bspl]/sigma;
y=sp(:)/sigma;

R = [zeros(size(Rb,1),Nm),Rb];
r=zeros(size(Rb,1),1);

D=diag([ones(Nm,1);zeros(Nb,1)]);
d=[zeros(Nm,1);-ones(Nb,1)];

if nargout>1
    [c_b,cov,dof,V]=contin(A,y,R,r,D,d,alpha);
    [tmp1,tmp2,dof0,V0]=contin(A,y,R,r,D,d,0);
    F = (V-V0)/V0*(length(sp)-dof0)/dof0;
    p=FTest(dof0,length(sp)-dof0,F);
    ec_b=sqrt(diag(cov));
    V=sqrt(V);
else
   % [tmp1,tmp2,tmp3,V]=contin(A,y,R,r,D,d,alpha);
    
     option=optimset('LargeScale','off','Display','off');
     [xi,V2]=lsqlin(A,y,D,d,[],[],[],[],[],option);
     V=sqrt(V2);
end