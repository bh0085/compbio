function z=my_and(x,y)

len=length(x);
z=x.*0-1;

for i=1:len
   if x(i)==1 & y(i)~=-1
       z(i)=1;
   elseif y(i)==1 & x(i)~=-1
       z(i)=1;
   elseif x(i)==0 & y(i)==0
       z(i)=0;
   end
    
    
end