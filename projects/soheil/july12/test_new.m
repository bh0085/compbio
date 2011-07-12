% test new

clc
clear all
close all
warning off

% control=0 % denovo network
% control=1  % sush net
%control=2  %  motif net
% control=3 %bind net
%control=4 % comb motif and bind

control=0;

a=['expression_c4d_n4_tt_24.mat'];
b=['results_test_24'];
c=['result_total_test_24'];


input_binary=struct('filename',a,'control',control,'results_filename',b,'iter_num',2,'f',1.1,'max_num3',5,'th_last',0.4,'output_name',c);

% a is the cluster name
% b is the results file name (not all the variables)
% c is a filename to save all variables
% iter_num shows the circuit complexity

binary_regression4(input_binary)