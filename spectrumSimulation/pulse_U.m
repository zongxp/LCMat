function U=pulse_U(s,pat,pw,B1max,phi,H0)
% U=pulse_U(s,pat,pw,B1max,phi[,H0])
% b as b1*gamma/2/pi
% angle from x axis
% phi: angle in the x-y plane
% H0: additional Hamiltonian applied during the pulse.
[rf,tp]=read_rf_varian(pat,pw,B1max);

if ~exist('H0','var')
    H0=0;
end

for i=1:length(rf)
    H=pulse_H(s,abs(rf(i)),angle(rf(i))*180/pi+phi)+H0;
    if i==1
     U=H2U(H,tp(i));
    else
     U=U*H2U(H,tp(i));
    end
end 


   
function [rf,tp]=read_rf_varian(rfname,pw,B1max)
%[rf,tp]=read_rf_varian(rfname,pw,B1max);

root=pwd;

    if isempty(rfname)
        a_tmp=repmat([0,0,1],[100,1]);     
    elseif strcmp(rfname,'x')
        a_tmp = [-pi/2,B1max,1];
    elseif strcmp(rfname,'y')
        a_tmp = [0,B1max,1];
    elseif strcmp(rfname,'-x')
        a_tmp = [pi/2,B1max,1];
    elseif strcmp(rfname,'-y')
        a_tmp = [pi,B1max,1];
    else
        a_tmp=textread(fullfile(root,rfname),'','commentstyle','shell');
        a_tmp(:,1)=a_tmp(:,1)*pi/180;
        a_tmp(:,2) = a_tmp(:,2)/max(a_tmp(:,2))*B1max;
    end
    
    
     t = a_tmp(:,3);
     tp = t/sum(t)*pw;
 rf = a_tmp(:,2).*(cos(a_tmp(:,1))+1i*sin(a_tmp(:,1)));