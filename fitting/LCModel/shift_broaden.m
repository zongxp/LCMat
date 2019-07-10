function y=shift_broaden(M,gm,sh,fstep)

 m=ifft(M);
 sw=fstep*length(M);
 t=(0:length(M)-1)/sw;
 Y=m.*exp(-(gm+1i*sh)*t');
 y=fft(Y);
 