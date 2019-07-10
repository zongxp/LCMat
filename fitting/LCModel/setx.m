function x0=setx(sh,s,gm)
% x0=setx(Ns,s,gm)
        
ngm=length(gm);
nsh=length(sh);
Ns=length(s)/2;

    if Ns>0
      x0(1:2*Ns)=s;
    end
    
    x0(2*Ns+1:2*Ns+ngm)=gm;
    ip=2*Ns+ngm;
    x0(ip+1:ip+nsh)=sh;
    