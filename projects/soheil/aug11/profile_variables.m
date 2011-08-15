function var=profile_variables(pf)

row_size=size(pf,1);

if row_size==1 & length(pf{1,1})==0
    var=[];
else
    var=[];
    for i=1:row_size
        var=union(pf{i,3},var); 
    end 
end