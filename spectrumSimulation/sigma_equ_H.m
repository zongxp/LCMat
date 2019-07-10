function sigma=sigma_equ_H(H)
% equilibrium density matrix for the zeeman energy
%H=spin_sys_H(s);
%H=H-s.f0*2*pi*spin_sys_U(s,'z'); 

[v,d]=eig(H);

        
U2=-d;
      
     % U2=(U2-trace(U2))/trace(abs(d));
     U2=U2-trace(U2);
     U2=U2/size(U2,1);
     U3=diag(diag(U2));
    
     
      sigma=v*U3*v';
      
      