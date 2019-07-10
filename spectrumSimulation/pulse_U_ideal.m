function U=pulse_U_ideal(s,fa,phi,ispin)
% U=pulse_U_ideal(s,fa,phi,ispin)
if exist('ispin','var')
 sx=spin_sys_U(s,'x',ispin);
 sy=spin_sys_U(s,'y',ispin);
else
 sx=spin_sys_U(s,'x');
 sy=spin_sys_U(s,'y');
end    

th=phi*pi/180; 
H=-(sx*cos(th)+sy*sin(th))*2*pi;

U=Uexp(H,1i*fa/360);