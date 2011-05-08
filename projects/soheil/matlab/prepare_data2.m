function [gene,tf,net_prob,gene_names,tf_names]=prepare_data2(filename)


% exp=open('expression_data.mat');
% exp=open('expression_c4d_n4_tt_4.mat');
% exp=open('expression_c4d_n4_tt_9.mat');
% exp=open('expression_c4d_n4_tt_18.mat');
exp=open(filename);
ntfs=open('network_tfs.mat');
wts=open('network_weights.mat');


gtfsnames=fieldnames(exp);
cexp = struct2cell(exp);

gene_names=fieldnames(ntfs); % the names are string
cntfs = struct2cell(ntfs);
cwts = struct2cell(wts);


g_num=size(gene_names,1);
gene=cell(g_num,1);

% gene i==gene_names(i)

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

net_prob=zeros(g_num,tf_num);

for i=1:g_num
    for j=1:size(cntfs{i},1)
        
        index=find_gene_index([cntfs{i,1}(j,:)],tf_names);
        net_prob(i,index)=cwts{i,1}(j,1);
    end
end


