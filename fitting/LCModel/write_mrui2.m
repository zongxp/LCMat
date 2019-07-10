function write_mrui2(data,prefix,sw)
% write_mrui2(data,prefix,sw)
% data is np*n_spectrum.
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
head(6)=127.74e6;
head(7)=3.0;
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





