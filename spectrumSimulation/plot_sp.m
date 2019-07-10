function ph=plot_sp(data,f_center,sw,phasecorr)
        
np=length(data);
x=f_center+linspace(-sw/2,sw/2,np);

fa=fft(data);
fa=fftshift(fa);
if phasecorr
 ph=phase(data(1))*180/pi;
else
 ph=0;
end
figure;plot(x,real(fa*exp(-1i*ph*pi/180)),'k-');
hold on;
plot(x,abs(fa),'b-');
%hold on;plot(x,imag(fa*exp(-1i*ph)),'r-');

set(gca,'XDir','reverse');
xlabel('Freq (ppm)');
%fac=fa*exp(-1i*ph);
%disp(sum(fac(2218:2405)));
%disp(sum(fac(1709:1843)));
%disp(sum(fac(1553:1708)));
