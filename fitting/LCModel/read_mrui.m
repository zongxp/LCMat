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


