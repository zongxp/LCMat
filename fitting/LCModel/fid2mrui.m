function fid2mrui(fid_prefix,ntrunc,ntrunc_end,fwhm,phase)
% fid2mrui(fid_prefix,ntrunc,ntrunc_end[,fwhm,phase])
% fid_prefix: cell array containing the names of the fid data folders
%  ntrunc: number of time points to remove at the beginning of the fid.
% ntrunc_end: number of time points to remove at the end of the fid
% fwhm: optional; additional broadening
% phase: optional; hard phase adjustment.

  fid_prefix=str2cell(fid_prefix);

  
for i=1:length(fid_prefix)  
  d=read_fid([fid_prefix{i},'.fid/fid']);
  sw=readPar(fid_prefix{i},'sw');
  d=d(ntrunc+1:end-ntrunc_end);
  d=conj(d);
 if exist('fwhm','var') && fwhm>0
    t=(0:length(d)-1)/sw;
    %tau = 1/fwhm/pi; %exponential decay  
    tau=sqrt(log(2))*4/fwhm/2/pi; %gaussian
    d=d.*exp(-(t'/tau).^2);
 end

 if exist('phase','var')
  d=d*exp(-1i*phase/180*pi);
 end  

  dall(:,i)=d;
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
function write_mrui2(data,prefix,sw)
% write_mrui2(data,prefix,sw)
% data is np*nt.
fid= fopen([prefix,'.mrui'],'w','ieee-be');
if ~exist('sw','var')
    sw=4000;
end
head=zeros(64,1);
head(1)=0;
head(2)=size(data,1);
head(3)=1/sw*1000;
head(4)=0;
head(5)=0;
head(6)=400e6;
head(7)=9.4;
head(8)=0;
head(9)=0;
head(10)=4.6;
head(14)=size(data,2);

fwrite(fid,head,'double');


d2=[real(data(:)),imag(data(:))];
d2=d2';
d2=d2(:);


fwrite(fid,d2,'double');

fclose(fid);


function z = read_fid(fid_dir,format)
% z = read_fid(fname[,format])
% the output is a complex matrix of dimension np/2,ntraces,nblocks.
% where np is the number of data points for each fid.
if exist([fid_dir,'.fid'],'dir') 
    fid_dir = [fid_dir,'.fid'];
 end

 if ~isdir(fid_dir) && exist(fid_dir,'file')
    fname=fid_dir;
 else
    
 if exist(fullfile(fid_dir,'fid.orig'),'file')
  fname=fullfile(fid_dir,'fid.orig');
 else
   fname=fullfile(fid_dir,'fid');
 end
 end


fid = fopen(fname,'r','ieee-be');

fh = read_rawfheader(fid);

z = zeros(fh.np/2,fh.ntraces,fh.nblocks);

status = dec2bin(fh.status);

for i=1:fh.nblocks
  
    for j=1:fh.nbheaders
        
     bh(i) = read_rawbheader(fid);
    end
    
    if status(end-3) == '1' && fh.ebytes == 4 && ~exist('format','var')
      format = 'float';
    elseif status(end-3) == '0' && fh.ebytes == 2 && ~exist('format','var') 
      format = 'int16';
    elseif status(end-3) == '0' && fh.ebytes == 4 && ~exist('format','var')
      format = 'int32';
    end
    
    tmp = fread(fid,fh.np*fh.ntraces,format);
    z(:,:,i) = reshape(tmp(1:2:end-1),fh.np/2,fh.ntraces)+1i*reshape(tmp(2:2:end),fh.np/2,fh.ntraces); 
    
   % dr(:,:,i) = reshape(tmp(1:end/2),fh.np/2,fh.ntraces); 
    %di(:,:,i) = reshape(tmp(end/2+1:end),fh.np/2,fh.ntraces); 
    
end
fclose(fid);


function fheader = read_rawfheader(fid)

f = fread(fid,6,'int32');
fheader.nblocks = f(1);
fheader.ntraces = f(2);
fheader.np = f(3);
fheader.ebytes = f(4);
fheader.tbytes = f(5);
fheader.bbytes = f(6);

a=fread(fid,2,'int16');

fheader.vers_id = a(1);
fheader.status = a(2);
fheader.nbheaders = fread(fid,1,'int32');


function val = readPar(fid_dir,field)
%val = readPar(fid_dir or parfile,field)



    
 if exist([fid_dir,'.fid'],'dir') 
    fid_dir = [fid_dir,'.fid'];
 elseif exist([fid_dir,'.img'],'dir') 
    fid_dir = [fid_dir,'.img'];
 end

 if ~isdir(fid_dir) && exist(fullfile(pwd,fid_dir),'file')
    fname=fid_dir;
 else
    
 if exist(fullfile(fid_dir,'procpar.orig'),'file')
  fname=fullfile(fid_dir,'procpar.orig');
 else
   fname=fullfile(fid_dir,'procpar');
 end
end

%a=textscan(fid,'%s%s%s%s%s%s%s%s%s%s%s','delimiter',' ');
%a=textread(parfile,'%s%s%s%s%s%s%s%s%s%s%s',);
a=textread(fname,'%s');

ind = strmatch(field,a,'exact');
if length(ind) ~=1
    val=[];
    warning('Field not found');
    return;
end

type = a{ind+1};

ind = ind+11;
na = str2num(a{ind});
if isempty(na)
    error('Number of values error');
end

val = a{ind+1};

         if ~isempty(str2num(val))
            val = zeros(1,na);
              for i=1:na
               val(i) = str2num(a{ind+i});
              end
              
              if type == '6'
                  val = val/1000000;
              end
         elseif na==1
             val = a{ind+1};
         else
             val=a(ind+1:ind+na)';
         end
            

