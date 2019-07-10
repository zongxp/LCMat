function y=myspline(j,n,t,k)
%  y=myspline(j,n,t,k)
% j: integer from 1 to length(k)-n-1
if j<1 || j>length(k)-n-1
    error('j out of range');
end

if n<0
    error('n should be a non-negative');
end

if n==0    
    y=zeros(size(t));
    for i=1:length(t)
     if t(i)>=k(j)&&t(i)<k(j+1)
      y(i)=1;
     end
    end
else    
    y=(t-k(j))/(k(j+n)-k(j)).*myspline(j,n-1,t,k);
    y=y+(k(j+n+1)-t)/(k(j+n+1)-k(j+1)).*myspline(j+1,n-1,t,k);
end

