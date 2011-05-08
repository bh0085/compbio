function m=mat_threshold(m,k)


 temp=sort(m(1:end),'descend');
 
 len=length(temp);
 max_val=[temp(1)];
 for i=2:len
     if temp(i)~=temp(i-1) & length(max_val)<k
         max_val=[max_val,temp(i)];
     end
 end

 
 for i=1:size(m,1)
     for j=1:size(m,2)
         if ismember(m(i,j),max_val)==0
            m(i,j)=0; 
         end
     end
 end
 
 
 m=m./sum(sum(m));