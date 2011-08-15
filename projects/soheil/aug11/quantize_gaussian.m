function y=quantize_gaussian(x)

m=mean(x);
sig=sqrt(var(x));
len=length(x);
y=x.*0;

for i=1:len
    if sig~=0
   y(i)=floor((x(i)-m)/sig+1/2); 
    else
    y(i)=0;    
    end
end