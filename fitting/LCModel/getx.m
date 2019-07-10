 function [sh,s,gm]=getx(x0,Ns,Nm,common,mac)
        
 
        
if common
    Nmp=1;
else
    if mac
        Nmp=Nm-1;
    else
        Nmp=Nm;
    end
end

    nx=2*Ns+Nmp+Nm;
   if length(x0)~=nx
       disp([Ns,Nmp,Nm,length(x0)]);
       
       error('Inconsistent input arguments');
   end
   
   if Ns>0
       s= x0(1:2*Ns);
   else
       s=[];
   end
   gm= x0(2*Ns+1:2*Ns+Nmp);
   ip=2*Ns+Nmp+1;
   sh= x0(ip:end);
  
   
   
   
   
   
   