function run_LCModel(params)


sbutton=get(gco,'String');

dname=get(params,'scan folder');
wdname=get(params,'alt. water folder');

extra_phase = get(params,'extra phase (deg)');

gain=readbPar(fullfile(dname,'acqp'),'RG');
nt = readbPar(fullfile(dname,'method'),'PVM_NAverages');  % no need to normalize by nt in PV6. smae for ntw.

TR = readbPar(fullfile(dname,'method'),'PVM_RepetitionTime');

if isempty(wdname)
    ws=readb_fid(fullfile(dname,'fid.refscan'));
    ntw=readbPar(fullfile(dname,'method'),'PVM_RefScanNA');
    gainw=readbPar(fullfile(dname,'method'),'PVM_RefScanRG');
    TRw= readbPar(fullfile(dname,'method'),'PVM_RepetitionTime');
     dsw=readbPar(fullfile(dname,'method'),'PVM_DigShift');
else
    ws=readb_fid(wdname);
    ntw=readbPar(fullfile(wdname,'method'),'PVM_NAverages'); % number of repetitions for acquiring water signal.
    gainw=readbPar(fullfile(wdname,'acqp'),'RG');  
    TRw= readbPar(fullfile(wdname,'method'),'PVM_RepetitionTime');
      dsw=readbPar(fullfile(dname,'method'),'PVM_DigShift');
end
np=size(ws(:),1);
ws=ws(:);
bl=mean(ws(end-round(np/10):end),1);
ws=ws(:)-repmat(bl,[np,1]);


TR = TR/1000;
TRw=TRw/1000;
fmrui='';%get(params,'mrui file');
if ~isempty(fmrui)
    [h,b]=read_mrui(fmrui);  %spectra preprocessed in mrui; 
else
    b=readb_fid(dname);
    ds=readbPar(fullfile(dname,'method'),'PVM_DigShift');
    
    b=b(ds+1:end);
    
end
b=b(:);
np=size(b,1);
bl=mean(b(end-round(np/10):end));
b=b(:)-bl;
 
[tmp,ind_max]=max(abs(ws));  % determine the first point for the water line.
a=b(:)/abs(mean(ws(ind_max:ind_max+3)));

fa=fft(a,[],1);
fa=fftshift(fa,1);

fa=fa*exp(-1i*extra_phase*pi/180);

if strcmp(get(gco,'String'),'Show Sp.')
    figure(101);plot(real(fa));
    return;
end

bl_ind=get(params,'baseline range');
bl_ind=str2num(bl_ind);

if strcmp(sbutton,'Calc Norm. Area')
    
    if length(bl_ind)~=4
        disp('Please define four points for the baseline (2 for left and right each)');
        return;
    end
   tmp =ts_detrend(real(fa),[bl_ind(1):bl_ind(2),bl_ind(3):bl_ind(4)],1);
   tmp=tmp/(length(ws)-ind_max+1)*ntw/nt*gainw/gain;  % area under the curve for the water spectrum is the first fid point times the number of data points for water line.
   
   figure;plot(bl_ind(1):bl_ind(4),tmp(bl_ind(1):bl_ind(4)),'b-');
   hold on;
   x = bl_ind([1,4]);
   plot(x,[0,0],'k-','LineWidth',0.5);
   x = bl_ind([1,2]);
   plot(x,[0,0],'r-','LineWidth',1);
   x = bl_ind([3,4]);
   plot(x,[0,0],'r-','LineWidth',1);
   
   res=sum(tmp(bl_ind(2)+1:bl_ind(3)-1));
   
   
   fprintf('Normalized peak area relative to water signal is %3.2e\n',res);    
   return;
end

name_met={};
if(get(params,'Cr')) name_met{end+1}='Cr'; end
if(get(params,'Glu')) name_met{end+1}='Glu'; end
if(get(params,'Gln')) name_met{end+1}='Gln'; end
if(get(params,'Lac')) name_met{end+1}='Lac'; end
if(get(params,'Ins')) name_met{end+1}='Ins'; end
if(get(params,'NAA')) name_met{end+1}='NAA'; end
if(get(params,'pCho')) name_met{end+1}='pCho'; end
if(get(params,'Tau')) name_met{end+1}='Tau'; end
if(get(params,'pCr')) name_met{end+1}='pCr'; end


dmod=get(params,'model sp. Folder');


M=get_model_spectra_local(name_met,dmod);

w=1;  %assume water content = 0.8 (unit is g/ml?)


gainr=gain/gainw;

sw=readbPar(fullfile(dname,'acqp'),'SW_h');
params=get_LCModel_params_local(0,3,w,nt,ntw,gainr,TR,TRw,name_met,sw);  %initialize parameters for LCModel.



fa_bl=real(fa(bl_ind(1):bl_ind(2)));
fa_bl=detrend(fa_bl,'linear');
s0=std(fa_bl);  % use flat part of the spectra to estimate the noise level.
 
params.std=s0;
params.last_is_mac=false;
[root,scanName]=fileparts(dname);

fname=sprintf('Results_%s',scanName);  % output file name.

params.alpha=[6,2]*3.8;
%params.alpha=[8,2]*5;

params.adjust_alpha=false;
%[10.5,2968];  %Initial regulaziation parameter.  These two parameters will be automatically adjusted in LCModel to achieve desired p value (0.1~<p~<0.8)
% and approximately equal regS (lineshape regularization term) and RegRb (baseline regularization term). 
% If automatic adjustment fail to reach termination condition, set
% params.adjust_alpha=false and adjust these two
% parameters manually.

LCModel(fa,params,M,fullfile(dname,fname));
  


function params=get_LCModel_params_local(tof,gm0,w,nt,ntw,dgain,TR,TRw,name_met,sw)
% params=get_params(tof,gm0,w,nt,ntw,dgain,TR)
% tof: shift of carrier frequency relative to water.
% gm0: initial linewidth broading parameter estimate (Units: Hz).
% w: water content divided by 0.8.
% nt: NEX for metabolite fid
% ntw: NEX for water fid
% dgain: gain difference between metabolite and water measurements.
% TR: repetition time.  Needed for accouting for the different steady-state magnetizations 
%     of water and metabolites (assumed T1 = 2 s and 1.4 s, respectively; values at 9.4 T. 
%     Modify accordingly at other field strength).
% 
if ~exist('dgain','var')
    dgain=30;
end
if ~exist('TR','var')
  TR=2.97;
end
FA=90;

f1=gre_transient(FA,nt,1,TR,1.4);
f2=gre_transient(FA,ntw,1,TRw,2);
f1=mean(f1);
f2=mean(f2);

params.name=name_met;
nspin_all=[5,6,5,5,4,6,6,13,4,5];
name_all={'Cr','GABA','Glu','Gln','Lac','Ins','NAA','pCho','Tau','pCr'};
c0_all=[10,3,10,4,1,3,10,5,10,1,450];

for i=1:length(name_met)
    
    ind=strmatch(name_met{i},name_all);
    nspin(i)=nspin_all(ind);
    c0(i)=c0_all(ind);
    
end

params.c0= c0; %initial concentration estimate.



%corr=f1/f2./(1000*110*(nt./nspin)/(w/0.8)*ntw/gain);

corr=ntw/nt*f2/f1*(1000*110)./nspin/dgain*(w/0.8);

params.corr=1./corr;

params.common=false;
params.fwhm=10;
gmi=gm0*ones(1,length(name_met));  % no extra broadening for Macromolecule
gmu=20*ones(1,length(name_met));
gml=zeros(1,length(name_met));
%gmi(6)=0;   %Myo-Inosital
%gmu(6)=6;
params.gamma=[gmi;gml;gmu];

params.tof=tof;   %shift from water;
params.frange=[1.0,4.25];

%sh0=tof+(4.675-3)*400;  %most model spectra except pCho and Ins are centered at 3 ppm.
%sh0=ones(1,10)*sh0;
%sh0(8)=tof+730;   %pCho model centered at 730 Hz upfield from water.
%sh0(6)=tof;    %Ins model centered at water

sh0=zeros(1,length(name_met));  %difference between experimental carrier frequence and that used for model spectra; assume they are the same

%sh0(end+1)=tof+730;  % Macromolecules

sh0(2,:)=sh0-20;
sh0(3,:)=sh0(1,:)+20;
params.sh0=sh0;     %difference between experimental carrier frequency and model spectrum carrier frequency. 
                   %to move the model spectrum to higher frequency to match
                   %the measured spectrum, sh0 should decrease and vice
                   %versa.


%range=1:2043;

params.f0=400e6;
params.sw=sw;

params.alpha=[44,723];  %1 for shape, 2 for baseline
params.adjust_alpha=true;
params.nofit=false;
params.noplot=false;

function s=gre_transient(b,N,M0,TR,T1)
% s=gre_transient(b,N,M0,TR,T1)
% N is the number of TRs.
% b is the flip angle in degree
% M0: initial magnetization


if ~exist('M0','var')
    M0=1;
end
s=zeros(1,N);
for i=1:N
  s(i)=M0*sin(b(1)*pi/180);
  M0=relax(M0*cos(b(1)*pi/180),TR,T1);
  
end

function s=relax(m0,delay,T1)

s=1-(1-m0)*exp(-delay/T1);



function M=get_model_spectra_local(a,dmod)
% [M,ra,offset]=get_model_spectra(a)
% offset frequency of the center of the model spectra relative to water in Hz.
if ~exist('a','var')
   a={'Cr','GABA','Glu','Gln','Lac','Ins','NAA','pCho','Tau','pCr'};
end

% for i=2:2%3:length(a)
%     
% special(a{i},0.004,0.02,4088,8013,4.675,4.675,true,0.001);
% end
% 

nm=length(a);

ra=zeros(1,nm);

for i=1:nm
    
    
  [hd,fid]=read_mrui(fullfile(dmod,a{i}));
 % fid=conj(fid);  
  fa=fft(fid);
  fa=fftshift(fa);
  ra(i)=abs(fid(1,:));
  ph=phase(fid(1));
  M(:,i)=fa*exp(-1i*ph);
end

M=M./repmat(ra,[size(M,1),1]);  %normalize to avoid numerical errors in error calculation for Macromolecules.
%ra=ones(1,nm);  % don't bother

function [head,s]=read_mrui(prefix)
%[head,s]=read_mrui(prefix)
% or s=read_mrui(prefix)

if length(prefix)<= 5 || ~strcmp(prefix(end-4:end),'.mrui')
  prefix=[prefix,'.mrui'];
end
fid=fopen(prefix,'r','ieee-be');
%fid=fopen([prefix,'.mrui'],'r','ieee-le');

a=fread(fid,64,'double');
N = a(2);
head.N = N;
head.freq = a(6);
head.ref_freqp=a(10);
head.ref_freqh=a(9);
head.dt = a(3);
head.t0=a(4);
head.phase0=a(5);
M=a(14);
head.M = a(14);

 b=fread(fid,2*N*M,'double');


s = b(1:2:end-1)+1i*b(2:2:end);


s=reshape(s,[N,M]);
fclose(fid);

if nargout==1
    head=s;
end

function [res,nzero]=readb_fid(d)


if ~exist(d,'dir')  % d is a file
    [d,fname,suf]=fileparts(d);
    fname=[fname,suf];
else
    fname='fid';
end

sz=readbPar(fullfile(d,'acqp'),'ACQ_size',true);
nr=readbPar(fullfile(d,'acqp'),'NR',true);
ni=readbPar(fullfile(d,'acqp'),'NI',true);

if numel(sz)>1
 ns = readbPar(fullfile(d,'acqp'),'NSLICES',true);
else
    ns=1;
end

fid=fopen(fullfile(d,fname),'r','ieee-le');

fmt=readbPar(fullfile(d,'acqp'),'GO_raw_data_format',false);

if strcmp(fmt,'GO_32BIT_SGN_INT')
 res=fread(fid,'int32');
 
else
    fprintf('%s :',fmt);
    fclose(fid);
    error('unknown format');
end

acqmod=readbPar(fullfile(d,'acqp'),'ACQ_experiment_mode',false);
if ~strcmp('SingleExperiment',acqmod)
    rcvrs=readbPar(fullfile(d,'acqp'),'ACQ_ReceiverSelect',false);
    nyes=strmatch('Yes',rcvrs);
    nch=length(nyes);    
   
else
    nch=1;
end

bs=readbPar(fullfile(d,'acqp'),'GO_block_size',false);

if strcmp('Standard_KBlock_Format',bs)
    sz1_old=sz(1)*nch;
 
 %s1=2.^ceil(log2(sz(1)*nch));
 %s2=2.^ceil(log2(sz(1)*nch/10))*10;
 %sz(1)=min(s1,s2)/nch;
 
 sz_ch=length(res)/(prod(sz(2:end))*ni*nr);
 %fprintf('Readout data points %d/%d\n',sz1_old,sz(1));
 nzero=sz_ch-sz1_old;
else
    sz1_old=sz(1)*nch;
    sz_ch=sz1_old;
    nzero=0;
    
end

nzero=nzero/2;
%res=res(1:2:end-1)+1i*res(2:2:end);
 
res=reshape(res,[2,sz_ch/2,sz(2:end)',ns,ni/ns*nr]);

res=res(1,:,:,:,:,:)+1i*res(2,:,:,:,:,:);
res=squeeze(res);
if size(res,1)~=1
 res=squeeze(res(1:sz1_old/2,:,:,:,:,:,:));
end
sz=size(res);
res=reshape(res,[sz(1)/nch,nch,sz(2:end)]);
fclose(fid);

function res=readbPar(fname,par,isnum)
%res=readbPar(fname,par,isnum)

if ~exist('isnum','var')
    isnum=true;
end

fid=fopen(fname,'r');
par=['##$',par,'='];
while 1
  b=fgetl(fid);
  if b==-1
      fclose(fid);
      error([par(4:end-1), ' not found']);
  end
  
  ind=strfind(b,par);
  if ~isempty(ind) 
    res=b(ind+length(par):end);
    break;
  end

end

if res(1)=='('
    sz=str2num(res(2:end-1));
     if numel(sz)==1
            sz=[1,sz];
     end
        
    if isnum 
     res=[];
     while 1
        
      b=fgetl(fid);
      if isnum
        tmp=str2num(b);
      else
       tmp=strread(b,'%s');  
      end
      res=[res,tmp];
      
      if length(res)==prod(sz) || ~isnum
          break;
      end
      
     end
    res=reshape(res,sz(end:-1:1));
    
    else
         b=fgetl(fid);
       res=strread(b,'%s');  
     
    end
else
    if isnum    
      res=str2double(res);
    end
    
    
end
fclose(fid);

function data = ts_detrend(data,tps,order)


   nt=length(data);
                    
                    ts = data(tps);
                    p = polyfit(tps,ts',order);
                    t = 1:nt;
                    
                    bl = zeros(1,nt);
                    for ip=0:order
                       bl = bl+t.^(order-ip)*p(ip+1);
                    end
                                           
                    data = data- bl'; 
            
        
