function bfid2mrui(fid_prefix,ntrunc,ntrunc_end,fwhm,phase)
% fid2mrui(fid_prefix,ntrunc,ntrunc_end[,fwhm,phase])

  fid_prefix=str2cell(fid_prefix);

  if ~exist('ntrunc','var')
      ntrunc=0;
  end
  if ~exist('ntrunc_end','var')
      ntrunc_end=0;
  end
  dall=[];
for i=1:length(fid_prefix)  
    
  nr=readbPar(fullfile(fid_prefix{i},'acqp'),'NR');
  d=readb_fid(fid_prefix{i});
  
  sw=readbPar(fullfile(fid_prefix{i},'acqp'),'SW_h');
  d=d(ntrunc+1:end-ntrunc_end,:);
 % d=conj(d);
 if exist('fwhm','var') && fwhm>0
    t=(0:size(d,1)-1)/sw;
    %tau = 1/fwhm/pi; %exponential decay  
    tau=sqrt(log(2))*4/fwhm/2/pi; %gaussian
    d=d.*repmat(exp(-(t'/tau).^2),[1,nr]);
 end

 if exist('phase','var')
  d=d*exp(-1i*phase/180*pi);
 end  

  dall(:,end+1:end+nr)=d;
end
 %d=d.*exp(-t'/tau);
 np=size(dall,1);
 bl=mean(dall(end-round(np/10):end,:),1);
 dall=dall-repmat(bl,[np,1]);
 
if length(fid_prefix)==1 
 write_mrui2(dall,fid_prefix{1},sw);
else
 write_mrui2(dall,[fid_prefix{1},'etc'],sw);
end    
    
%%

