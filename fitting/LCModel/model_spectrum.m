function [y_base,M2]=model_spectrum(s,M,gm,sh,fstep,mac)
% y_base=model_spectrum(s,M,gm,sh,fstep)
% s: 1*(2*Ns+1)
% M: np*9, model spectra
% gm: 1*9, broadening linewidth (FWHM) in Hertz;
% sh: 1*9, shift in Hertz.
% spectra width in Hertz.
Nm=size(M,2);
if length(gm)==1
    gm=gm*ones(1,Nm);
end
y_base=zeros(size(M));

Ns=(length(s)-1)/2;
M2=zeros(size(M));
for i=1:Nm
    if i~=Nm || ~mac   
      M2(:,i)=shift_broaden(M(:,i),gm(i)*pi,sh(i)*2*pi,fstep);
      M3=conv(M2(:,i),s);
      y_base(:,i)=M3(Ns+1:end-Ns);
    else
      y_base(:,i)=shift_broaden(M(:,i),0,sh(i)*2*pi,fstep);
    end

    
end
