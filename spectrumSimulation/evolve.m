function out=evolve(sigma,H,tp)

%[v,d]=eig(H);
%U2=v*exp(1i*d*tp)*v';
U2=Uexp(H,tp*1i);
out=U2'*sigma*U2;