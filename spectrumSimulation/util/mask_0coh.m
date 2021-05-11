function mask=mask_0coh(s)

n=s.nspins;
mask=false(2^n,2^n);
for i=1:2^n
  for j=1:2^n
      si=dec2bin(i-1,n);
      sj=dec2bin(j-1,n);
      
      if length(find(si=='0')) == length(find(sj=='0'))
          mask(i,j)=1;
          
      end
  end
end

  