% itercluster analysis
clc
clear all
close all

warning off

addpath(genpath('C:\My Doc\Projects\comp bio\Matlab-Network framework\to send\network modules-lassoMCMC\glmnet_matlab'))
addpath(genpath('C:\My Doc\Projects\comp bio\Matlab-Network framework\to send\network modules-lassoMCMC\data'))
addpath(genpath('C:\My Doc\Projects\comp bio\Matlab-Network framework\to send\network modules-lassoMCMC\network data'))
addpath(genpath('C:\My Doc\Projects\comp bio\Matlab-Network framework\to send\network modules-lassoMCMC\data2'))


a=['expression_c4d_n4_intercluster.mat'];
b=['results_intercluster.mat'];
c=['results_total_intercluster.mat'];

% control=0 % denovo network
% control=1  % sush net
% control=2  %  motif net
% control=3 %bind net
% control=4 % comb motif and bind

control=4;

input_binary=struct('filename',a,'control',control,'results_filename',b,'iter_num',2,'f',1.1,'max_num3',10,'th_last',0.4,'output_name',c);
% input_binary=struct('filename',a,'results_filename','results','iter_num',1,'f',1.1,'max_num3',5,'output_name','result_binary_12');
binary_regression4(input_binary)

 load(b);
 load(c);
        
 infered_network;
 quality_test;
 
 
 figure
 hist(quality_test)
 
 rep_cell=module_repetition(network_cell);
 
 
[sorted_module_rep,I]=sort([rep_cell{:,2}],'descend');

index=1;


while rep_cell{I(index),2}>=3
    index=index+1;
end

index

if index>=2
rep_cell2=rep_cell(I(1:index-1),:);
else
    rep_cell2=rep_cell(I(1),:);
end
% rep_cell{I(1:2),1}

figure
bar(sorted_module_rep)
title('module repetition-intercluster')

figure
hist([rep_cell2{:,2}])
title('hist of module repetition-intercluster')

save('rep_cell_inter','rep_cell2')

