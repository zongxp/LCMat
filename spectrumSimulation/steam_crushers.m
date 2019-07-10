function steam_crushers(mname,te,tm,T2s,np,sw,freq)
% steam(mname,te,tm,T2s,np,sw,freq)
% mname: name of the matlab function that defines the spin system.
% te: echo time.  Units: s
% tm: time between the second and the third rf pulses. Units: s
% T2s: decay time constant for the observed signal.  Note FWHM is 1/T2s/pi.
% np: number of complex data points to acquire
% sw: spectral width
% freq: Frequency of the rf pulses and receiver (units: ppm)
% 
%

fname = sprintf('SpinOperators_%s.mat',mname);
if ~exist(fname,'file')

    fprintf('Calculating Hamiltonian and spin operators for %s\n',mname);
    s=eval(mname);
    sigeq=sigma_equ(s);
    [H0,HJ,Hz]=spin_sys_H(s);
   
    s.sx=spin_sys_U(s,'x');
    s.sy=spin_sys_U(s,'y');
    s.sz=spin_sys_U(s,'z');
    save(fname);
    fprintf('Hamiltonian and spin operators for %s saved in %s\n',mname,fname);
else
   load(fname,'s','H0','sigeq');
   fprintf('Using Hamiltonian and spin operators for %s saved in %s\n',mname,fname);
end


%t1=[0 2 0 2 0 2 0 2,0 2 0 2 0 2 0 2,1 3 1 3 1 3 1 3,1 3 1 3 1 3 1 3];
%t2=[0 0 0 0 2 2 2 2,0 0 0 0 2 2 2 2,1 1 1 1 3 3 3 3,1 1 1 1 3 3 3 3];
%t3=[0 0 2 2 0 0 2 2,0 0 2 2 0 0 2 2,1 1 3 3 1 1 3 3,0 0 2 2 0 0 2 2];
%t7=[0 2 2 0 2 0 0 2,0 2 2 0 2 0 0 2,1 3 3 1 3 1 1 3,1 3 3 1 3 1 1 3];

dy=s.sy;
dx=s.sx;
%f=linspace(0,500*(1-1/nt),nt);
%f2=linspace(0,500*(1-1/nt2),nt2);
%f=linspace(0,2000*(1-1/nt),nt);
%f2=linspace(0,1000*(1-1/nt2),nt2);
f=[0,2500];
f2=[0,2500,5000,7500];  %four needed to cancel double-quantum coherences.
nt=length(f);
nt2=length(f2);

data=zeros(np,nt*nt2);

U=pulse_U_ideal(s,90,90);

for i=1:nt%length(t1)
  for j=1:nt2
   disp([f(i),f2(j)]);
   Hsp=-f(i)*2*pi*s.sz;
   
   sig=evolveU(sigeq,U);
   
   if 1
   sig=evolve(sig,H0,te/2-0.0001);
   sig=evolve(sig,H0+Hsp,0.0001);
   %U=pulse_U_ideal(s,90,t2*90);
  
   sig=evolveU(sig,U);
 %  sig=evolve(sig,HJ+Hz,0.030);
 %  sig=crusherU(sig,s.sz);    %the results using crusher gradients are
   %close to phase cycling, but not exactly the same (why ?).probably Hz is
   %not correctly calculated.
   
  % sig=evolve(sig,H0,tm-0.0001);
  % Hsp2=-f2(j)*2*pi*s.sz;
  % sig=evolve(sig,H0 + Hsp2,0.0001);
    sig=evolve(sig,H0,tm);
    Hsp2=-f2(j)*2*pi*s.sz;
    sig=evolve(sig,Hsp2,0.0001);   %[H0,Hsp2]=0
   
   
   sig=evolveU(sig,U);
   sig=evolve(sig,H0,te/2-0.0001);
   sig=evolve(sig,H0+Hsp,0.0001);
   end
   
   detect=(dx-1i*dy);
   data_t=fid(s,sig,H0,detect,freq*s.f0/1e6,sw,np);
   t=(0:np-1)/sw;
   afwhm=1/T2s/pi;
   data_t=data_t.*exp(-t*afwhm*pi);
   data(:,i+(j-1)*nt)=data_t;
  end
%subplot(4,6,i);
%plot_sp(data_t,3,0);

end
save(mname);
fprintf('Results saved in %s.mat\n',mname);
figure;
ph=plot_sp(mean(data,2),freq,sw*1e6/s.f0,1);
write_mrui2(-mean(data*exp(1i*ph),2),mname,sw);
    
fprintf('Simulated spectrum saved in %s.mrui\n',mname); 
    
      