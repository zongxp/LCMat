 function y=Rs(n)
        
        y=zeros(n+2,n);
        
        for i=1:n
            y(i,i)=1;
            y(i+1,i)=-2;
            y(i+2,i)=1;
        end
        