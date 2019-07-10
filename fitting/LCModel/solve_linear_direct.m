function [V,c_b,ec_b,pval]=solve_linear_direct(sp,model_base,bspl,Rb,alpha,sigma)
%[c_b,ec_b,p,V]=solve_linear_direct(sp,model_base,bspl,Rb,alpha,sigma)

Nm=size(model_base,2);
Nb=size(bspl,2);
A=[model_base,bspl]/sigma;
y=sp(:)/sigma;
X = [A;[zeros(size(Rb,1),Nm),alpha*Rb]];
y=[y;zeros(size(Rb,1),1)];

D=diag([ones(Nm,1);zeros(Nb,1)]);
d=[zeros(Nm,1);ones(Nb,1)];
%option=optimset('LargeScale','on','Display','off');
option = optimoptions('lsqlin','Display','off','Algorithm','active-set');
 
c_b=lsqlin(X,y,-D,d,[],[],[],[],[],option);
V = y-X*c_b;

if nargout > 2
    
    binding=find(c_b(1:Nm)<1e-6);
    X(:,binding)=[];
    
    covy=diag([ones(1,length(sp)),zeros(1,size(Rb,1))]);
    Xinv=inv(X'*X)*X';
    cov=Xinv*covy*Xinv';
    
    ec_b=zeros(size(c_b));
    inz=setdiff(1:length(c_b),binding);
    ec_b(inz)=sqrt(diag(cov));
    dof0=Nm+Nb;
    
    X0 = [A;zeros(size(Rb,1),Nm+Nb)];
    c_b0=lsqlin(X0,y,-D,d,[],[],[],[],[],option);
   % c_b0=lsqlin(X0,y,-D,d,option);
    V0 = y-X0*c_b0;
    Vnorm=sum(V.^2);
    V0norm=sum(V0.^2);
    F = (Vnorm-V0norm)/V0norm*(length(sp)-dof0)/dof0;
    pval=FTest(dof0,length(sp)-dof0,F);
    
end


function p = FTest(dof1, dof2, F, dof2max)
%
% p = FTest(dof1, dof2, F, <dof2max>)
%
% Computes p-value given F-value. p and F can be vectors.
% dof1 = dof of numerator (Number of Rows in RM)
% dof2 = dof of denominator (DOF)
%
% Ref: Numerical Rec in C, pg 229.
%
%
%


%
% FTest.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:29 $
%    $Revision: 1.3 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
%



if(nargin ~= 3 & nargin ~= 4)
  msg = 'Usage: p = FTest(dof1, dof2, F, <dof2max>)';
  qoe(msg);error(msg);
  return;
end

if(length(dof1) > 1)
  error('dof1 must be a scalar');
end
if(length(dof2) > 1)
  error('dof2 must be a scalar');
end

if(exist('dof2max') ~= 1) dof2max = []; end
if(isempty(dof2max))      dof2max = dof2; end
dof2 = min(dof2,dof2max);

if F<0
    p=0;
return;
end

z = dof2./(dof2 + dof1 * F);
p = betainc(z, dof2/2, dof1/2);

return;
