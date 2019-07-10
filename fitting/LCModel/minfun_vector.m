function y=minfun_vector(x,sp,M,bspl,Ns,Nb,common,fstep,alpha,sigma,gm0)
% x should be a row vector.
% alpha: regulizaor 1*2
% sigma: errors: 1*3.  First for spectrum. 2nd: linewidth. 3rd: shift.
Nm=size(M,2);
np=length(sp);

[c,sh,b,s,gm]=getx(x,Ns,Nb,Nm,common);

if Ns>0
    s(end+1)=1-sum(s);
else
    s=1;
end

model=model_spectrum(c,s,M,gm,sh,fstep);

bl=baseline(bspl,b);

t1=real(sp-bl-model)/sigma(1);
t2=(Rs(Ns*2+1)*s')*alpha(1);
t3=Rs(np)*bl;
t3=t3(3:end-2)*alpha(2);
t4=(gm-gm0)/sigma(2);
t5=sh/sigma(3);

if Ns==0
    t2=[];
end

y=[t1(:);t2(:);t3(:);t4(:);t5(:)];

   