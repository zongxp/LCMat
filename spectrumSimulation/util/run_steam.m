function run_steam(params)


 a={'Cr','GABA','Glu','Gln','Lac','Ins','NAA','pCho','Tau','pCr'};
 
 doutput=get(params,'output folder');
 
 scand=get(params,'scan folder');
 f0=get(params,'Freq (MHz)');
 
 te = readbPar(fullfile(scand,'method'),'PVM_EchoTime');
 tm=readbPar(fullfile(scand,'method'),'StTM');
 T2s = 0.05;
 np=readbPar(fullfile(scand,'method'),'PVM_DigNp');
 shift=readbPar(fullfile(scand,'method'),'PVM_DigShift');
  
sw=readbPar(fullfile(scand,'acqp'),'SW_h');
 
freq_rf=4.675+readbPar(fullfile(scand,'method'),'PVM_SpecOffsetHz')/f0;
freq_rcv=freq_rf;

 mkdir(doutput);
 cd(doutput);
 
 ideal_pulse=true;
 pw=0.0027;  % change tomorrow;
 
 for i=1:length(a)
   steam(a{i},te/1000,tm/1000,T2s,np-shift,sw,freq_rf,freq_rcv,ideal_pulse,pw,f0)

 end
 
 cd ..;

