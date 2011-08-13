%test

clc
clear all
close all
warning off


a= '/home/unix/benh/compbio/data/soheil/expression_c4d_n4_tt_65.mat'
b= '/home/unix/benh/compbio/data/soheil/mat_out_0/expression_c4d_n4_tt_65.mat'
c= '/home/unix/benh/compbio/data/soheil/mat_all_out_0/expression_c4d_n4_tt_65.mat'


% control=0 % denovo network
% control=1  % modencode net
% control=2  %  motif net
% control=3 %bind net
% control=4 % comb motif and bind

control=0;

%iter_num: is the determining the maximum circuit complexity
% prior: penalty on complexity
% a,b,c: input output file names
input_binary=struct('filename',a,'control',control,'results_filename',b,'iter_num',2,'f',1.1,'max_num3',10,'th_last',0.2,'output_name',c,'prior',0.75);


binary_regression5(input_binary);


