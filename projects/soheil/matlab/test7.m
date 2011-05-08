% test 7
clc
clear all
close all
warning off

input=struct('filename','expression_c4d_n4_tt_4.mat','out_iter_num',2,'in_iter_num',5,'k',6,'beta',4,'f_mix',1,'f_sim',0.9,'f_en_in',0.95,'f_en_out',0.95,'th_cor',0.5,'trunc_value',3,'degree_bound',3,'output_name','result');

%out_iter_num=50
%in_iter_num=1000

%k=4:10
% beta=4:k beta is less or equal to k

%f_mix=0.5:0.1:2
%f_sim= 0.5:0.05:0.95
%f_en_in=0.5:0.05:1
%f_en_out=0.5:0.05:1
%th_cor=0.4:0.05:0.8
%trunc_value=2:5
%degree_bound=3:6

MCMC_noonlinear(input);
