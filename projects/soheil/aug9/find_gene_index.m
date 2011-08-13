function index=find_gene_index(name,list)

% list is a cell of size (g_num+tf_num)*1
% name is a char array

n=size(list,1);
stop=0;
if n>1
    for i=1:n
        if stop==0
            if sum((list{i,1}==name)-1)==0
                index=i;
                stop=1;
            else
                index=0;
            end
        end
    end
else
    n=size(list,2);
    for i=1:n
        if stop==0
            if sum((list{1,i}==name)-1)==0
                index=i;
                stop=1;
            else
                index=0;
            end
        end
    end
end

