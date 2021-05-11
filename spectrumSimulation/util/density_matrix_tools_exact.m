function density_matrix_tools_exact()

mname='GABA';
if 0
s=Glu();
%s=proton;
%s=threeproton();
%[rf,tp]=read_rf_varian('P10.RF',1e-3,250/0.15);

sigeq=sigma_equ(s);
[H0,HJ,Hz]=spin_sys_H(s);
   
s.sx=spin_sys_U(s,'x');
s.sy=spin_sys_U(s,'y');
s.sz=spin_sys_U(s,'z');

save(['density_matrix_tools_',mname]);
end

load(['density_matrix_tools_',mname]);

np=8198;
dy=s.sy;
dx=s.sx;
t1=0;
t7=0;
U=pulse_U_ideal(s,90,t1*90);
x=[1,1.111];
for i=1:2
   sig=evolveU(sigeq,U);
   if 1  %change to 0 for fid
   sig=evolve(sig,H0,0.002);
  % sig=crusherU(sig,s.sz);
   %U=pulse_U_ideal(s,90,t2*90);
   
   sig=evolveU(sig,U);
 %  sig=crusherU(sig,s.sz);
   sig=evolve(sig,H0,0.03);
   sig=evolve(sig,s.sz*2000,x(i));
   %U=pulse_U_ideal(s,90,t3*90);
   sig=evolveU(sig,U);
   sig=evolve(sig,H0,0.002);
   %sig=crusherU(sig,s.sz);
   
   end
   
   if i==1
       sig1=sig;
   else
       sig2=sig;
   end
end


   sw=8000;
   detect=(dx-1i*dy)*exp(-1i*t7*90);
   data_t=fid(s,sig,H0,detect,1200,sw,np);
   t=(0:np-1)/sw;
   afwhm=8;
   data_t=data_t.*exp(-t*afwhm*pi);


figure;
plot_sp(data_t,3,1,true);


%write_mrui2(data_t(1:2043),mname);
    

    
      