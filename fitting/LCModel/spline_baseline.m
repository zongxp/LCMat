function s=spline_baseline(t,k,p)
% k: knots
% p: control points. length(k)-4
m=length(k);
n=3;  %cubic spline
s=zeros(size(t));
for i=1:m-n-1
    s=s+p(i)*myspline(i,n,t,k);
end




