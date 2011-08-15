function [gene,tf,net_prob,gene_names,tf_names]=prepare_data5(filename,control)

% control=0 % denovo network
% control=1  % sush net
%control=2  %  motif net
% control=3 %bind net
%control=4 % comb motif and bind

% determining gene_names and tf_names

ntfs=open('sushnet.mat');
gene_names=fieldnames(ntfs); % the names are string
cntfs = struct2cell(ntfs);

tf_names=cell(1,1);
for i=1:size(cntfs,1)
    for j=1:size(cntfs{i,1},1)
        if i==1 & j==1
            tf_names{1,1}= [cntfs{1,1}(1,:)];
        elseif find_gene_index([cntfs{i,1}(j,:)],tf_names)==0
            tf_names=[tf_names,[cntfs{i,1}(j,:)]];
        end
    end
end
tf_names=tf_names';

% loading expression levels

exp=open(filename);
gtfsnames=fieldnames(exp);
cexp = struct2cell(exp);

g_num=size(gene_names,1);
gene=cell(g_num,1);

for i=1:g_num
    index=find_gene_index(gene_names{i,1},gtfsnames);
    gene{i,1}=cexp{index,1}';
end

tf_num=size(tf_names,1);
tf=cell(tf_num,1);

% gene i==gene_names(i)

for i=1:tf_num
    index=find_gene_index(tf_names{i,1},gtfsnames);
    tf{i,1}=cexp{index,1}';
end


check_names=ones(g_num,tf_num);
% zero in this matrix means they are the same genes
for i=1:g_num
    for j=1:tf_num
        temp_names=gene_names{i,:}==tf_names{j,:};
        if ismember(0,temp_names)==0
            check_names(i,j)=0;
        end
    end
end


%********
%loading the network

net_prob=zeros(g_num,tf_num);

if control==0
    for i=1:g_num
        for j=1:tf_num
            if check_names(i,j)==1
                net_prob(i,j)=1;
            end
        end
    end
    
elseif control==1
    network='sushnet.mat';
    for i=1:g_num
        for j=1:tf_num
            
            gene1=gene_names{i,1};
            gene2=tf_names{j,1};
            
            e=network_edge(network,gene1,gene2);
            if e==1
                net_prob(i,j)=e;
            end
        end
    end
elseif control==2
    network='mnet.mat';
    for i=1:g_num
        for j=1:tf_num
            
            gene1=gene_names{i,1};
            gene2=tf_names{j,1};
            
            e=network_edge(network,gene1,gene2);
            if e==1
                net_prob(i,j)=e;
            end
        end
    end
    
elseif control==3
    network='bnet.mat';
    for i=1:g_num
        for j=1:tf_num
            
            gene1=gene_names{i,1};
            gene2=tf_names{j,1};
            
            e=network_edge(network,gene1,gene2);
            if e==1
                net_prob(i,j)=e;
            end
        end
    end
    
elseif control==4
    
    network='mnet.mat';
    for i=1:g_num
        i
        for j=1:tf_num
            
            gene1=gene_names{i,1};
            gene2=tf_names{j,1};
            
            e=network_edge(network,gene1,gene2);
            if e==1
                net_prob(i,j)=e;
            end
        end
    end
    
    network='bnet.mat';
    for i=1:g_num
        i
        for j=1:tf_num
            
            gene1=gene_names{i,1};
            gene2=tf_names{j,1};
            
            e=network_edge(network,gene1,gene2);
            if e==1
                net_prob(i,j)=e;
            end
        end
    end
    
end



