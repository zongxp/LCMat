function H=pulse_H(s,b,offset,phi)
% H=pulse_H(s,b,offset,phi)
% b as b1*gamma/2/pi
% angle from x axis
% phi: angle in the x-y plane
% offset: rf field offset
% returned hamiltonian is in the rotating frame at 0 ppm.

sx=s.sx;
sy=s.sy;
sz=s.sz;
th=phi*pi/180;

H=-(sx*cos(th)+sy*sin(th))*b*2*pi-sz*offset*2*pi;


   