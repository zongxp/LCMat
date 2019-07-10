function steam(mname,te,tm,T2s,np,sw,freq_rf,freq_rcv,ideal_pulse,pw,f0)
% steam(mname,te,tm,T2s,np,sw,freq_rf,freq_rcv,ideal_pulse,pw)
% mname: name of the matlab function that defines the chemical shifts and J couplings of the metabolite.
%        Several metabolites have been defined in the spectrum_simulation
%        directory.
% te: echo time.  Units: s
% tm: time between the second and the third rf pulses. Units: s
% T2s: Apparent T2 (units: s).  Note FWHM is 1/T2s/pi.
% np: number of complex data points.
% sw: spectral width in Hz.
% freq_rf: Carrier Frequency of the rf pulses. (units: ppm. water is 4.675 ppm.)
% freq_rcv: Reference frequency of the receiver. (default: 4.675 ppm);  
% ideal_pulse: whether to use ideal hard pulse or the shaped "P10" and
% "P11" pulses.
% pw: pulse width of the shaped pulse (Units: s). Dimension: 1*1 or 1*3.
% If 1*1, it will be set to pw =[pw*1.26,pw,pw];
% f0: resonance freq in MHz.
% A program for simulating metabolite spectra in a STEAM sequence.
% Author: Xiaopeng Zong; NeuroImaging Lab; University of Pittsburgh.
% Last update 5/27/2013.

root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'spectrum_simulation'));
if ~exist('freq_rcv','var')  || isempty(freq_rcv)
    freq_rcv=4.675;
end

fname = sprintf('SpinOperators_%s.mat',mname);
if ~exist(fname,'file')  %check whether spin Hamiltonian already calculated before.

    fprintf('Calculating Hamiltonian and spin operators for %s\n',mname);
    s=eval(sprintf('%s(%f)',mname,f0*1e6));
   
    [H0,HJ,Hz]=spin_sys_H(s);
    sigeq=sigma_equ(s);
    s.sx=spin_sys_U(s,'x');
    s.sy=spin_sys_U(s,'y');
    s.sz=spin_sys_U(s,'z');
    save(fname);
    fprintf('Hamiltonian and spin operators for %s saved in %s\n',mname,fname);
else
   load(fname,'s','H0','sigeq');
   fprintf('Using Hamiltonian and spin operators for %s saved in %s\n',mname,fname);
end


dy=s.sy;
dx=s.sx;

f=[0,2500];  %two resonant freq due to crusher gradients in the TE period.

nt=length(f);
f2=0;  %Only need to calculate for one frequency during TM period.  Directly set single and higher order coherences to 0.
nt2=length(f2);
data=zeros(np,nt*nt2);
Mz=zeros(1,nt*nt2);

H0_rf=H0+s.sz*freq_rf*s.f0*1e-6*2*pi;  %the saved H0 are on the rotating reference frame at 0 ppm.
H0_rcv=H0+s.sz*freq_rcv*s.f0*1e-6*2*pi; 
if ideal_pulse
  U1=pulse_U_ideal(s,90,90);
  U2=pulse_U_ideal(s,90,180);
  U3=pulse_U_ideal(s,90,270);
  d1=0;
  d2=0;
  d3=0;
else
if length(pw)==1
  pw =[pw*1.26,pw,pw];
end
asym=0.17;
U1=pulse_U(s,'spectrum_simulation/P10.RF',pw(1),1.67/pw(1),90,H0_rf);
U2=pulse_U(s,'spectrum_simulation/P11.RF',pw(2),1.67/pw(2),180,H0_rf);
U3=pulse_U(s,'spectrum_simulation/P10.RF',pw(3),1.67/pw(3),270,H0_rf);

d1=(pw(1)+pw(2))*asym;
d2=(pw(2)+pw(3))*(1-asym);
d3=pw(3)*asym;

end

for i=1:nt%length(t1)
  for j=1:nt2
   Hsp=-f(i)*2*pi*s.sz;
% 1st pulse
   sig=evolveU(sigeq,U1);
% evolve
   sig=evolve(sig,H0_rf,te/2-d1);
   sig=evolve(sig,Hsp,0.0001);
% 2nd pulse  
   sig=evolveU(sig,U2);
% evolve
   sig(~mask_0coh(s))=0;
   sig=evolve(sig,H0_rf,tm-d2);
% 3rd pulse  
   sig=evolveU(sig,U3);
% evolve   
   sig=evolve(sig,Hsp,0.0001);
   sig=evolve(sig,H0_rf,te/2-d3);
% acquire data   
   detect=(dx-1i*dy);
   data_t=fid(s,sig,H0_rcv,detect,0,sw,np);
   t=(0:np-1)/sw;
   afwhm=1/T2s/pi;
   data_t=data_t.*exp(-t*afwhm*pi);
   data(:,i+(j-1)*nt)=data_t;
   Mz(i+(j-1)*nt)=trace(sig*s.sz);
  end

end

save(mname);
fprintf('Results saved in %s.mat\n',mname);

ph=plot_sp(mean(data,2),freq_rcv,sw*1e6/s.f0,1);
write_mrui2(-mean(data*exp(1i*ph*pi/180),2),mname,sw);
    
fprintf('Simulated spectrum saved in %s.mrui\n',mname); 

function write_mrui2(data,prefix,sw)
% write_mrui2(data,prefix,sw)
% data is np*nt.
fid= fopen([prefix,'.mrui'],'w','ieee-be');
if ~exist('sw','var')
    sw=4000;
end
head=zeros(64,1);
head(1)=0;
head(2)=size(data,1);
head(3)=1/sw*1000;
head(4)=0;
head(5)=0;
head(6)=400e6;
head(7)=9.4;
head(8)=0;
head(9)=0;
head(10)=4.6;
head(14)=size(data,2);

fwrite(fid,head,'double');


d2=[real(data(:)),imag(data(:))];
d2=d2';
d2=d2(:);


fwrite(fid,d2,'double');

fclose(fid);

    
      