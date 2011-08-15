% network data
clc
clear all 
close all

my_path='C:\My Doc\Projects\comp bio\Matlab-Network framework\to send\network modules-lassoMCMC\network data';
addpath 'C:\My Doc\Projects\comp bio\Matlab-Network framework\to send\network modules-lassoMCMC\network data';

% net=open('bnet.mat');
% all_gene_list=fieldnames(net);
% net_cell = struct2cell(net);

% index1=find_gene_index(gene_1,all_gene_list);
% index2=find_gene_index(gene_2,all_gene_list);

gene1='FBgn0027563';
gene2='FBgn0013263';

network='C:\My Doc\Projects\comp bio\Matlab-Network framework\to send\network modules-lassoMCMC\network data\bnet.mat';
e=network_edge(network,gene1,gene2)
tf_list=network_children(network,gene1)