function data=fid(s,sigma,H0,detect,offset,sw,np)
     
      sz=s.sz;
      data=zeros(1,np);
      H0=H0+offset*2*pi*sz;
      U=H2U(H0,1/sw);
     for i=1:np
         data(i)=trace(sigma*detect);
         sigma=evolveU(sigma,U);
     end