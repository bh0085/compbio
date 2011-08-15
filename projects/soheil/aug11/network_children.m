function tf_list_cell=network_children(network,gene_1)

% this function takes a network, if there is an edge between gene_2 to
% gene_2, e=1 o.w. e=0

net=open(network);
all_gene_list=fieldnames(net);
net_cell = struct2cell(net);

index1=find_gene_index(gene_1,all_gene_list);
index1;

if index1>=1
    tf_list=net_cell{index1}.tfs;
    m=size(tf_list,1);
    tf_list_cell=cell(m,1);
    for i=1:m
        
       tf_list_cell{i,1}=tf_list(i,:); 
    end
   
    
else
    tf_list_cell=cell(1,1); % error, gene_1 is not in the list
end

