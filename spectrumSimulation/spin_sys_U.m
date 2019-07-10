function U=spin_sys_U(s,xyz,ind_spin)
%  U=spin_sys_U(s,xyz[,ind_spin])
% default sum of all spins in ind_spin.
% spin operator of a multiple spin system.

n=s.nspins;
U=zeros(2^n,2^n);
so=spin_U(1/2,xyz);

if ~exist('ind_spin','var')
    ind_spin=1:n;
end
% first spin on the most significant bit.

for i=1:length(ind_spin)
    
Utmp=zeros(2^n,2^n);
  for j=1:2^n
    for k=1:2^n
        
    str1=dec2bin(j-1,n);
    str2=dec2bin(k-1,n);
    
    state1=str2double(str1(ind_spin(i)));
    state2=str2double(str2(ind_spin(i)));
     
     if state1==0
        lv=[1,0];
     else
        lv=[0,1];
     end
     if state2==0
        rv=[1;0];
     else
        rv=[0;1];
     end
    
     same_for_others=1;
     for l=1:n
       state1=str2double(str1(l));
       state2=str2double(str2(l));
       if l~=ind_spin(i) && state1~=state2  
        same_for_others=0;
        break;
       end
     end
     
     Utmp(j,k)=lv*so*rv*same_for_others;
    end
  end
  
  U=Utmp+U;
end

