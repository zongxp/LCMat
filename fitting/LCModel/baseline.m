function [bl,bl_base]=baseline(bspl,b)
% b is the control point for spline

bl=zeros(size(bspl,1),1);
bl_base=zeros(size(bspl,1),length(b));

for i=1:length(b)
    bl=bl+b(i)*bspl(:,i);
    bl_base(:,i)=bspl(:,i);
end

