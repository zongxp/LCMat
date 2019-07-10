function [dydgm,dydsh]=shift_broaden_derivative(M,gm,sh,fstep)
% derivative of gm and sh
 m=ifft(M);
 sw=fstep*length(M);
 t=(0:length(M)-1)/sw;
 %Y=m.*exp(-(gm+1i*sh)*t');
 %y=fft(Y);
 
 dYdgm=-t'.*m.*exp(-(gm+1i*sh)*t');
 dydgm=fft(dYdgm);
 
 dydsh=1i*dydgm;
 