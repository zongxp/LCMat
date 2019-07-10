Iz=spin_sys_U(Lac,'z');
[H,H_lab]=spin_sys_H(Lac);

sig=zeros(size(Iz));
sig(1,3)=1;


t=linspace(0,2*pi);

sig2=zeros([size(Iz),100]);
for i=1:length(t)
    sig2(:,:,i)=Uexp(-Iz,t(i)*1i)*sig*Uexp(Iz,t(i)*1i);
end

figure;plot(squeeze(real(sig2(1,3,:))));

