function [xfit,pval,regS,regB]=LCModel(sp,p,M,outf)
%[xfit,chi2,p_chi2]=LCModel(sp,p,M,outf)
% sp: measured spectra N*1
% M: model spectrum  N*m.  m is the number of metabolites.
% p: parameter structure containing:
%   range: indices array to select the range in sp and M for fitting.
%   corr: coversion factor.  sp = sum(M.*corr.*c0,2); 
%   name: name of the model spectra.  1*m string cells
%   sh0: 3*m, initial, lower and upper limits of metabolite frequency shift in rad/s.
%   outf: mat file name to save output
%   sw: spectral width.
%   f0: resonance frequency
%   last_is_mac: the last model spectrum is macrmolecule, no broadening
%   will be applied on this one.  default: true.
%   tof: carrier frequency shift relative to water. Hz.
% outf: name of the output file.
% Author: Xiaopeng Zong (zongxp@gmail.com), NeuroImaging Laboratory, University of Pittsburgh.
% Last updated: 5/28/2013.

sh0=p.sh0;
corr=p.corr;
name = p.name;
%f0=400e6;
f0=p.f0;
sw=p.sw;
fwhm=p.fwhm;
frange=p.frange;
if isfield(p,'last_is_mac')
    mac=p.last_is_mac;
else
    mac=true;
end

tof=p.tof;
np=length(sp);
if np~=size(M,1)
    error('Data point mismatch between model and measured spectra');
end

freq=((1:np)-np/2-1)*sw/np+tof;
freq=freq/400+4.675;
indl= find(freq(1:end-1)<frange(1) & freq(2:end)>frange(1));
indh= find(freq(1:end-1)<frange(2) & freq(2:end)>frange(2));

range=indl:indh;
freq=freq(range);
%sw=8000;
% 0.25 ppm minimum spacing between knots.
npfit=length(range);
fstep=sw/length(sp);
dk=max([fwhm*1.5,50]);
Nb=round(fstep*npfit/dk);

Ns=round(1.5*fwhm/fstep/2);  %line shape convolution
%Ns=0;
s=ones(1,2*Ns+1);
if Ns>0
 for i=1:2*Ns+1
  s(i)=exp(-(i-Ns-1).^2/Ns*2*log(2));
 end
end
s=s/sum(s);

Nm=size(M,2);
sigma=zeros(1,3);
sigma(1)=p.std;   % error for the spectrum
sigma(2)=4;    %broadening 
sigma(3)=0.004*f0*1e-6;

if isfield(p,'common')
    common=p.common;
else
    common = false;   %use common broadening for all metabolites (gm)
end


gmm=p.gamma;
gm0=gmm(1,:);  %extra broadening
if common  && size(gmm,2)~=1
      error('broadening factor should have only one column');
else
    if mac && size(gmm,2)~=Nm-1
       error('Broadening factor should have %d columns',Nm-1);
    elseif ~mac && size(gmm,2)~=Nm
       error('Broadening factor should have %d columns',Nm);
    end
end
              
x0=setx(sh0(1,:),s(1:end-1),gmm(1,:));
xl=setx(sh0(2,:),zeros(1,2*Ns),gmm(2,:));
xu=setx(sh0(3,:),ones(1,2*Ns),gmm(3,:));

alpha=p.alpha;

%xfit=fmincon(@(x) minfun(x,sp(range),M(range,:),Ns,Nb,common,sw,alpha,sigma,gm0),x0,A,B,Aeq,Beq);
nk=Nb+4;
dk=(npfit-1)/(nk-2);
k=linspace(1-dk,npfit+dk,nk);
bspl=zeros(npfit,Nb);

for i=1:Nb
 bspl(:,i)=myspline(i,3,1:npfit,k);
end
tic;
 R1=Rs(npfit);
 Rb=Householder(R1(3:end-2,:)*bspl);
%Rb=R1(3:end-2,:)*bspl;


pvala=[];
a2array = [];  %old alpha2 values 
aratio = [];  %ratio between alpha(2) and alpha(1)
rRatio = [];  %ratio between regB/regS
while 1
if ~p.nofit

   %  [tmp,g]=minfun_r2(x0,sp(range),M(range,:),bspl,Rb,Ns,common,fstep,alpha,sigma,gm0,sh0(1,:));
   %  Jstr=double(g~=0);
     options=optimset('MaxIter',50,'Display','iter','UseParallel','always','TolFun',0.0001,'Jacobian','off');%,'JacobPattern',Jstr);
    % xfit=lsqnonlin(@(x2) minfun_r2(x2,sp(range),M(range,:),bspl,Rb,Ns,common,fstep,alpha,sigma,gm0),x0,xl,xu,options);
    xfit=lsqnonlin(@(x2) minfun_r2(x2,sp(range),M(range,:),bspl,Rb,Ns,common,fstep,alpha,sigma,gm0,sh0(1,:),mac),x0,xl,xu,options);
else
   xfit=x0;
end  

fprintf('\n fitting time = %f s\n\n',toc);

%calculate regS and regB
   [sh,s,gm]=getx(xfit,Ns,Nm,common,mac);  %x only contain non-linear parameters.
   vfinal= minfun_r2(xfit,sp(range),M(range,:),bspl,Rb,Ns,common,fstep,alpha,sigma,gm0,sh0(1,:),mac);
   nRb=size(Rb,1);     
   regB=sum(vfinal(npfit+1:npfit+nRb).^2)/Nb;
  if Ns==0
       regS=regB;
  else
       regS =sum(vfinal(npfit+nRb+1:npfit+nRb+2*Ns+3).^2)/2/Ns; 
  end

% calculate errors
   if Ns>0
     s(end+1)=1-sum(s);
   else
     s=1;
   end
   model_base=model_spectrum(s,M,gm,sh,fstep,mac);

[V,c_b,ec_b,pval]=solve_linear_direct(real(sp(range)),real(model_base(range,:)),bspl,Rb,alpha(2),sigma(1));
  fprintf('pval = %3.2f,regS = %4.3f regRb = %4.3f, alpha =[%4.3f,%4.3f] \n\n',pval,regS,regB,alpha(1),alpha(2));

  if p.nofit
      c_b(1:Nm)=p.c0(1,:).*corr;
      break;
  end

if (regB > regS/1.5 && regB < regS*1.5)&& (pval<=0.9 && pval>=0.1)
  break;
else
    if ~p.adjust_alpha
      break;
    else
     [alpha,a2array,aratio,pvala,rRatio]=set_alpha(regB,regS,pval,alpha,a2array,aratio,pvala,rRatio);
    end
end
    
      fprintf('Try with alpha =[%4.3f,%4.3f] \n',alpha(1),alpha(2));

end

bl=bspl*c_b(Nm+1:end);
spfit=model_base(range,:)*c_b(1:Nm)+bl;


   fprintf('name       conc(mM)           shift(Hz) broadening (Hz)\n');
  if common
      gm=gm*ones(1,Nm-1);
  end
  if mac
    Nmp=Nm-1;
  else
    Nmp=Nm;
  end

   for i=1:Nmp
      fprintf('%4s        %4.1f (%4.1f)      %4.1f       %4.1f\n',...
      name{i},c_b(i)/corr(i),ec_b(i)/corr(i),sh(i),gm(i));
   end
   
   if mac
      fprintf('%4s        %4.1f (%4.1f)       %4.1f \n',...
      name{Nm},c_b(Nm)/corr(Nm),ec_b(Nm)/corr(Nm),sh(Nm));
   end
   
%% plot results.
%load('LCModel');

save(outf);

if ~p.noplot
    
    if Ns>0
     figure;
     s(end+1) = 1-sum(s);
     ns=length(s);
     plot((-ns/2+1:ns/2)*fstep,s);
     xlabel('Freq (Hz)');
     title('Line shape Kernal');
    end

yl=[-0.3,1.1]*max(real(sp(range)));
%{
figure;
plot(freq,real(sp(range)-model_base(range,end)*c_b(Nm)-bl),'k');
set(gca,'XDir','reverse');
xlim([freq(1),freq(end)]);
ylim(yl);
title('spectrum after subtracting out baseline');
%}
figure;
subplot(3,1,1);
plot(freq,real(sp(range)),'k');
hold on;
plot(freq,real(spfit),'r');
plot(freq,-0.1*max(real(sp(range)))+real(sp(range)-spfit),'b');
set(gca,'XDir','reverse');

xlim([freq(1),freq(end)]);
ylim(yl);
legend('Data','Fit','Residual');
title('Spectrum');
subplot(3,1,2);

%name = {'Cr','GABA','Glu','Gln','Lac','Ins','NAA','PCho','Tau'};
symb={'k-','b-','r-','b-','k-','k-','k-','b-','k-','c-','g-'};


for i=1:Nmp    % the last one for macromolecules is displayed in the third plot
  plot(freq,real(model_base(range,i))*c_b(i),symb{i});
  hold on;
end

set(gca,'XDir','reverse');
xlim([freq(1),freq(end)]);
ylim(yl);
title('Individual Metabolites');
subplot(3,1,3);
plot(freq,real(bl),'b-');
hold on;
if mac
  plot(freq,real(model_base(range,end))*c_b(Nm),'g-');
  plot(freq,real(model_base(range,end))*c_b(Nm)+real(bl),'k-');
  legend('spline','Mac','Total');
else
    legend('spline');
end
title('Baseline Fit');
xlabel('Freq. (ppm)');
%hold on;
%plot(real(bl),'k');
%plot(real(bl+sp_Mac),'k');
set(gca,'XDir','reverse');

xlim([freq(1),freq(end)]);
ylim(yl);
set(gcf,'Position',[1         127         560         671]);
saveas(gcf,[outf,'.fig']);
end

