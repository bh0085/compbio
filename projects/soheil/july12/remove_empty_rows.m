function pf_new=remove_empty_rows(pf)

row_size=size(pf,1);

if length([pf{:,3}])==0
    pf_new=cell(1,5);
else
    
    pf_new=cell(1,5);
    
     if length([pf{1,3}])~=0
      pf_new=[pf_new;pf(1,:)]; 
   end
    
   for i=2:row_size-1 
    if length([pf{i,3}])==0
    pf_new=[pf(1:i-1,:);pf(i+1:end,:)];
        
    end
   end
   
   if length([pf{end,3}])~=0
      pf_new=[pf_new;pf(end,:)]; 
   end
    

end


if size(pf_new,1)>=2
    
   pf_new=pf_new(2:end,:); 
end


