function z=my_not(x)

len=length(x);
z=x.*0;

for i=1:len
   if x(i)==1 
       z(i)=-1;
   elseif x(i)==-1
       z(i)=1;
   end
    
end