function press()

mname='Lac';
if 1
s=Lac();
%s=proton;
%s=threeproton();
%[rf,tp]=read_rf_varian('P10.RF',1e-3,250/0.15);

sigeq=sigma_equ(s);
%H0=spin_sys_H(s);
   [H0,HJ,Hz]=spin_sys_H(s);
   
s.sx=spin_sys_U(s,'x');
s.sy=spin_sys_U(s,'y');
s.sz=spin_sys_U(s,'z');

save(['density_matrix_tools_',mname]);
end

load(['density_matrix_tools_',mname]);

np=4000;
t1=0;
t2=0;
t3=0;
t7=0;
%t1=[0 2 0 2 0 2 0 2,0 2 0 2 0 2 0 2,1 3 1 3 1 3 1 3,1 3 1 3 1 3 1 3];
%t2=[0 0 0 0 2 2 2 2,0 0 0 0 2 2 2 2,1 1 1 1 3 3 3 3,1 1 1 1 3 3 3 3];
%t3=[0 0 2 2 0 0 2 2,0 0 2 2 0 0 2 2,1 1 3 3 1 1 3 3,0 0 2 2 0 0 2 2];
%t7=[0 2 2 0 2 0 0 2,0 2 2 0 2 0 0 2,1 3 3 1 3 1 1 3,1 3 3 1 3 1 1 3];

dy=s.sy;
dx=s.sx;
%nt=2;
%nt2=2;
%f=linspace(0,500*(1-1/nt),nt);
%f2=linspace(0,500*(1-1/nt2),nt2);
%f=linspace(0,2000*(1-1/nt),nt);
%f2=linspace(0,1000*(1-1/nt2),nt2);
%f=[0,250];
f=0;
%f2=[0,500];
f2=0;  %spoiling gradient is not neccessary for simulating PRESS sequence.
te1=linspace(0.002,0.14,25);
te2=0.004;
nt=length(f);
nt2=length(f2);
data=zeros(np,nt*nt2);

U=pulse_U_ideal(s,90,t1*90);

U2=pulse_U_ideal(s,180,t1*90);

for i=1:length(te1)   
   sig=evolveU(sigeq,U);
   sig=evolve(sig,H0,te1(i)/2);
   
   sig=evolveU(sig,U2);
   sig=evolve(sig,H0,te1(i)/2+te2/2);
   
   sig=evolveU(sig,U2);
   sig=evolve(sig,H0,te2/2);
   
   
   sw=6000;
   detect=(dx-1i*dy)*exp(-1i*t7*90);
   data_t=fid(s,sig,H0,detect,1200,sw,np);
   t=(0:np-1)/sw;
   afwhm=8;
   data_t=data_t.*exp(-t*afwhm*pi);
   data(:,i)=data_t;
end

save(mname);

ang=angle(data);
figure;plot(ang(1,:));

figure;

%mname='Asp';
%load(mname);
 %t=(0:np-1)/sw;
 %  afwhm=3;
   
%data=data.*repmat(exp(t'*afwhm*pi),[1,4]);
ph=plot_sp(mean(data,2),3,1);
fprintf('phase = %f\n',ph);
%write_mrui2(-mean(data(1:2043,:)*exp(1i*ph),2),[mname,'_5Hz']);
    
%write_mrui2(-mean(data(1:2043,:)*exp(1i*ph),2),mname);
 
    
      