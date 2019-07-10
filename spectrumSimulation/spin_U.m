function U=spin_U(s,xyz)
% spin operators

if s==1/2
 if strcmp(xyz,'x')
    U=1/2*[0,i;-i,0];
 elseif strcmp(xyz,'y')
    U=1/2*[0,1;1,0];
 elseif strcmp(xyz,'z')
    U=1/2*[1,0;0,-1];
 end
end