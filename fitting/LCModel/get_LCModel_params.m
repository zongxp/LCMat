function params=get_LCModel_params(tof,gm0,w,nt,ntw,dgain,TR)
% params=get_params(tof,gm0,w,nt,ntw,dgain,TR)
% tof: shift of carrier frequency relative to water.
% gm0: initial linewidth broading parameter estimate (Units: Hz).
% w: water content divided by 0.8.
% nt: NEX for metabolite fid
% ntw: NEX for water fid
% dgain: gain difference between metabolite and water measurements in db.
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
f1=Mss(FA,TR,1.4);
f2=Mss(FA,TR,2);


[M,ra,offset]=get_model_spectra_local();
params.name={'Cr','GABA','Glu','Gln','Lac','Ins','NAA','PCho','Tau','pCr'};
nspin=[5,6,5,5,4,6,6,13,4,5];
gain=10.^(dgain/20);
corr=f1/f2./(1000*110*(ra/nt./nspin)/(w/0.8)*ntw/gain);
params.corr=corr;

params.common=false;
params.fwhm=10;
gmi=gm0*ones(1,10);  % no extra broadening for Macromolecule
gmu=10*ones(1,10);
gml=zeros(1,10);
gmi(6)=0;   %Myo-Inosital
gmu(6)=6;
params.gamma=[gmi;gml;gmu];

params.tof=tof;   %shift from water;
params.frange=[1.0,4.25];

%sh0=tof+(4.675-3)*400;  %most model spectra except pCho and Ins are centered at 3 ppm.
%sh0=ones(1,10)*sh0;
%sh0(8)=tof+730;   %pCho model centered at 730 Hz upfield from water.
%sh0(6)=tof;    %Ins model centered at water

sh0=tof-offset;

%sh0(end+1)=tof+730;  % Macromolecules

sh0(2,:)=sh0-20;
sh0(3,:)=sh0(1,:)+20;
params.sh0=sh0;     %difference between experimental carrier frequency and model spectrum carrier frequency. 
                   %to move the model spectrum to higher frequency to match
                   %the measured spectrum, sh0 should decrease and vice
                   %versa.


%range=1:2043;
params.c0=[10,3,10,4,1,3,10,5,10,1,450];  %initial concentration estimate.

params.f0=400e6;
params.sw=4000;

params.alpha=[44,723];  %1 for shape, 2 for baseline
params.adjust_alpha=true;
params.nofit=false;
params.noplot=false;


function s=Mss(fa,TR,T1)
%s=ss_GEEPI(fa,TR[,T1,TE,T2s])

fa=fa*pi/180;
e1=exp(-TR/T1);
s=sin(fa)*(1-e1)/(1-cos(fa)*e1);

