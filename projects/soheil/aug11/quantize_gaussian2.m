function y=quantize_gaussian2(x,th)

m=mean(x);
sig=sqrt(var(x));
len=length(x);
y=x.*0;

for i=1:len
    if sig~=0
        if (x(i)-m)/sig>th
            y(i)=1;
        elseif (x(i)-m)/sig<-th
            y(i)=-1;
        else
            y(i)=0;
        end
    else
        y(i)=0;
    end
    
    
end