function [gene,tf,net_prob,gene_names,tf_names]=prepare_data6(filename,control)

% control=0 % denovo network
% control=1  % sush net
%control=2  %  motif net
% control=3 %bind net
%control=4 % comb motif and bind

% determining gene_names and tf_names

a=open('gene_names.mat');
gene_names=a.gene_names;

b=open('tf_names.mat');
tf_names=b.tf_names;


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



%********
%loading the network

if control==0
    net=open('denovonet_bdtnp.mat');
    net=net.net_prob;
elseif control==1
    net=open('sush_bdtnp.mat');
    net=net.net_prob;
elseif control==2
    net=open('mnet_bdtnp.mat');
    net=net.net_prob;
elseif control==3
    net=open('bnet_bdtnp.mat');
    net=net.net_prob;
elseif control==4
    net=open('combnet_bdtnp.mat');
    net=net.net_prob;
end

net_prob=net;
    
    
    
    
    
    
    
    
    
    
    