function gmrui_callback(params)

fname=get(params,'fid file');

[dname,fid_name,suf]=fileparts(fname);

fid_name=[fid_name,suf];

ntrunc=readbPar(fullfile(dname,'method'),'PVM_DigShift');
ntrunc_end=0;


% fid2mrui(fid_prefix,ntrunc,ntrunc_end[,fwhm,phase])

  d=readb_fid(fname);
  
  sw=readbPar(fullfile(dname,'acqp'),'SW_h');
  d=squeeze(d(ntrunc+1:end-ntrunc_end));

   d=d(:);

  
 %d=d.*exp(-t'/tau);
 np=size(d,1);
 bl=mean(d(end-round(np/10):end,:),1);
 d=d-repmat(bl,[np,1]);
 
 write_mrui2(d,fullfile(dname,fid_name),sw);

    
%%

