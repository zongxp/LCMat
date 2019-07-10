function [y,g]=minfun_r2(x,sp,M,bspl,Rb,Ns,common,fstep,alpha,sigma,gm0,sh0,mac)
% y=minfun_r2(x,sp,M,bspl,Rb,Ns,common,fstep,alpha,sigma,gm0)
% x: nonlinear parameters to be determined. should be a row vector.
% sp: data
% M: model spectra
% bspl: baseline splines
% Rb: regularizor for baseline.  Should be a square matrix of size(bspl,2)
% Ns: number of points for lineshape function
% common: whether use a common broadening for all metabolites
% alpha: regulizor 1*2
% fstep: frequency increment in the spectrum
% sigma: errors: 1*3.  First for spectrum. 2nd: linewidth. 3rd: shift.
% gm0: expected R2
% mac: whether the last model spectrum is macromolecules.  No broadening
% will be applied for this spectrum.
Nm=size(M,2);
np=length(sp);
[sh,s,gm]=getx(x,Ns,Nm,common,mac);  %x only contain non-linear parameters.

if Ns>0
    s(end+1)=1-sum(s);
    %s=s./sum(s);
else
    s=1;
end

[model_base,M2]=model_spectrum(s,M,gm,sh,fstep,mac);
%V=solve_linear(real(sp),real(model_base),bspl,Rb,alpha(2),sigma(1));

[V,c_b]=solve_linear_direct(real(sp),real(model_base),bspl,Rb,alpha(2),sigma(1));
% Rb: length(sp)-2*Nb
if Ns==0
    t2=[];
else
    t2=(Rs(Ns*2+1)*s')*alpha(1);
end

t4=(gm-gm0)/sigma(2);
t5=(sh-sh0)/sigma(3);

y=[V(:);t2(:);t4(:);t5(:)];

% calculate Jacobian
if nargout == 1 
    return;
end

if common
    Nmp=1;
else
    Nmp=Nm;
end
nx=2*Ns+Nmp+Nm;
g=zeros(length(y),nx);
%
if Ns>0
  Rsm=  Rs(Ns*2+1);
  g(np+1:np+size(Rsm,1),1:size(Rsm,2))=alpha(1)*Rsm;
end

for l=1:Nmp
  g(np+length(t2)+l,2*Ns+l) = 1/sigma(2); 
end


for l=1:Nm
  g(np+length(t2)+Nmp+l,2*Ns+Nmp+l) = 1/sigma(3); 
end

PHI=zeros(np,length(c_b));
PHI(:,1:Nm)=real(model_base);
PHI(:,Nm+1:end) = bspl;


dPHI=zeros(np,length(c_b),nx);
% derivative wrt to S
if Ns>0
    k=Ns:np-Ns-1;
    for n=1:2*Ns
        dPHI(k,1:Nm,n)=M2(k-n+1+Ns,:);
    end
end

% derivative wrt gamma and shift
for i=1:Nm 
    if common
        [Mtmp,Mtmp2] = shift_broaden_derivative(M(:,i),gm*pi,sh(i)*2*pi,fstep);
        M3=conv(Mtmp,s);
        dPHI(:,i,2*Ns+1)=M3(Ns+1:end-Ns)*pi;
        M4=conv(Mtmp2,s);
        dPHI(:,i,2*Ns+1+i)=M4(Ns+1:end-Ns)*2*pi; 
    else       
        [Mtmp,Mtmp2]= shift_broaden_derivative(M(:,i),gm(i)*pi,sh(i)*2*pi,fstep);
        M3=conv(Mtmp,s);
        dPHI(:,i,2*Ns+i)=M3(Ns+1:end-Ns)*pi;
        M4=conv(Mtmp2,s);
        dPHI(:,i,2*Ns+Nm+i)=M4(Ns+1:end-Ns)*2*pi;
    end
end

for i=1:nx
    dph=real(dPHI(:,:,i));
    g(1:np,i) = -(eye(np)-PHI*(inv(PHI'*PHI)*PHI'))*dph*c_b; 
end







