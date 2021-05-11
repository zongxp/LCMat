function [H,H_lab,HJ,Hz]=spin_sys_H(s,cf)

% cf: carrier frequency in ppm. Default 0
% Hamiltonian in the rotating frame with frequency cf (units: ppm)

if ~exist('cf','var')
    cf=0;
end

n=s.nspins;
Ja=s.J;
pair=s.pairs;
offset = s.shift;
f0=s.f0;

offset=(offset-cf)*f0*1e-6;

if size(pair,1)~=length(Ja)  || length(offset)~=n
    error('spin system error');
end

H=zeros(2^n,2^n);
HJ=zeros(2^n,2^n);

for i=1:n
    sz=spin_sys_U(s,'z',i);
    H=H-sz*offset(i)*2*pi;
end
szall=spin_sys_U(s,'z');
Hz=H+szall*mean(offset)*2*pi;

for i=1:length(Ja)
      p1=pair(i,1);
      p2=pair(i,2);
      sx1=spin_sys_U(s,'x',p1);
      sx2=spin_sys_U(s,'x',p2);
      sy1=spin_sys_U(s,'y',p1);
      sy2=spin_sys_U(s,'y',p2);
      sz1=spin_sys_U(s,'z',p1);
      sz2=spin_sys_U(s,'z',p2);
      HJ=HJ+Ja(i)*(sx1*sx2+sy1*sy2+sz1*sz2)*2*pi;   
end

H=H+HJ;

H_lab=H-szall*f0*2*pi;

