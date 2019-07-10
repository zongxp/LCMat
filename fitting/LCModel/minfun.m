function y=minfun(x,sp,M,Ns,Nb,common,fstep,alpha,sigma,gm0)
% x should be a row vector.
% alpha: regulizaor 1*2
% sigma: errors: 1*3.  First for spectrum. 2nd: linewidth. 3rd: shift.
Nm=size(M,2);
np=length(sp);


[c,sh,b,s,gm]=getx(x,Ns,Nb,Nm,common);

model=model_spectrum(c,s,M,gm,sh,fstep);
bl=baseline(np,b);
t1=sum(real(sp-bl-model).^2)/sigma(1)^2;
t2=sum((Rs(Ns*2+1)*s').^2)*alpha(1)^2;
t3=Rs(np)*bl;
t3=sum(t3(3:end-3).^2)*alpha(2)^2;
t4=sum((gm-gm0).^2)/sigma(2)^2;
t5=sum(sh.^2)/sigma(3)^2;

y=t1+t2+t3+t4+t5;
