function [alpha,a2array,aratio,pvala,rRatio]=set_alpha(regB,regS,pval,alpha,a2array,aratio,pvala,rRatio)

pvala(end+1)=pval;
aratio(end+1)=alpha(2)/alpha(1);
a2array(end+1)=alpha(2);


rRatio(end+1)=regB/regS;

if length(find(~isnan(rRatio)))>=2  && any(rRatio>1) && any(rRatio<1)  
   rtmp = interp1(rRatio(~isnan(rRatio)),aratio(~isnan(rRatio)),1);
else
    rtmp = aratio(end)/rRatio(end);
end
  
if length(pvala)>=2 && any(pvala>0.5) && any(pvala<0.5) 
    alpha(2)=interp1(pvala,a2array,0.5);
elseif pval>0.5
    alpha(2)=alpha(2)*(1+0.5*(pval-0.5)/0.5);
else
    alpha(2)=alpha(2)/(1+0.5*(0.5-pval)/0.5);
end


if pval>=0.05
  alpha(1)=alpha(2)/rtmp;
else
    alpha(1)=alpha(1)/(1+0.5*(0.5-pval)/0.5);
    rRatio(end)=NaN;
end



function yi=my_interp(x,y,xi)


[x,i]=sort(x);
y=y(i);

if xi<=x(1)
    X = x(1:2);
    Y=y(1:2);
elseif xi>=x(end)
    X = x(end-1:end);
    Y=y(end-1:end);
else
    ind=find(xi>=x(1:end-1) & xi<=x(2:end));
    X = x([ind(1),ind(1)+1]);
    Y = y([ind(1),ind(1)+1]);
end

X=[X(:),ones(length(X),1)];
kb=inv(X'*X)*X'*Y(:);
yi=kb(1)*xi+kb(2);




