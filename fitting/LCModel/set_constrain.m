function [A,B,Aeq,Beq]=set_constrain(Nm,Ns,Nb,common)
        
% dimension of x: (Nm)+(2*Ns+1)+(1|Nm)+Nm+Nb        

if common
    nx=(Nm)+(2*Ns+1)+1+Nm+Nb;
    A=zeros(Nm+1,nx);
    B=zeros(Nm+1,1);
    Aeq=zeros(1,nx);
    Beq=1;
    for i=1:Nm
        A(i,i)=-1;
    end
    A(Nm+1,2*Ns+Nm+2)=-1;
    Aeq(Nm+1:2*Ns+Nm+1)=1;
else
    nx=(Nm)+(2*Ns+1)+Nm+Nm+Nb;
    A=zeros(Nm*2,nx);
    B=zeros(Nm*2,1);
    Aeq=zeros(1,nx);
    Beq=1;
    for i=1:Nm
        A(i,i)=-1;
        A(Nm+i,2*Ns+Nm+1+i)=-1;
    end
    Aeq(Nm+1:2*Ns+Nm+1)=1;
end
