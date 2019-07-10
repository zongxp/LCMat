   function U3=Uexp(U,factor)
  % calculate exp(U*factor)
  % U needs to be a Hermitian matrix
  % this is not always equal to Uexp(U*factor,1);
      [v,d]=eig(U);
        
      U2=exp(d*factor);
      
      U3=zeros(size(U));
      for i=1:size(U3,1)
          U3(i,i)=U2(i,i);
      end
      
      U3=v*U3*v';