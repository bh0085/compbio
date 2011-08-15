function [gene,tf,net_prob,gene_names,tf_names]=prepare_data4(filename,control)

% control=0 % denovo network
% control=1  % sush net
%control=2  %  motif net
% control=3 %bind net
%control=4 % comb motif and bind

exp=open(filename);
ntfs=open('sushnet.mat');
gtfsnames=fieldnames(exp);
cexp = struct2cell(exp);
gene_names=fieldnames(ntfs); % the names are string
cntfs = struct2cell(ntfs);

g_num=size(gene_names,1);
gene=cell(g_num,1);



for i=1:g_num
    index=find_gene_index(gene_names{i,1},gtfsnames);
    gene{i,1}=cexp{index,1}';
end

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
    
    
    ntfs2=open('sushnet.mat');
    cntfs2 = struct2cell(ntfs2);
    
    
    for i=1:g_num
        for j=1:size(cntfs2{i,1},1)
            
            index=find_gene_index([cntfs2{i,1}(j,:)],tf_names);
            if index~=0
                net_prob(i,index)=1;
            end
        end
    end
    
    
elseif   control==2 | control==3
    
    if control==2
        ntfs2=open('mnet.mat');
        cntfs2 = struct2cell(ntfs2);
        
    elseif control==3
        
        ntfs2=open('bnet.mat');
        cntfs2 = struct2cell(ntfs2);
    end
    
    
    for i=1:g_num
        for j=1:size(cntfs2{i,1}.tfs,1)
            
            index=find_gene_index([cntfs2{i,1}.tfs(j,:)],tf_names);
            if index~=0
                net_prob(i,index)=1;
            end
        end
    end
    
elseif control==4
    
    ntfs2=open('mnet.mat');
    cntfs2 = struct2cell(ntfs2);
    
    for i=1:g_num
        for j=1:size(cntfs2{i,1}.tfs,1)
            
            index=find_gene_index([cntfs2{i,1}.tfs(j,:)],tf_names);
            if index~=0
                net_prob(i,index)=1;
            end
        end
    end
    
    ntfs3=open('bnet.mat');
    cntfs3 = struct2cell(ntfs3);
    
    for i=1:g_num
        for j=1:size(cntfs3{i,1}.tfs,1)
            
            index=find_gene_index([cntfs3{i,1}.tfs(j,:)],tf_names);
            if index~=0
                net_prob(i,index)=1;
            end
        end
    end
    
end




