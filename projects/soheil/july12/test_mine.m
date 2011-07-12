function input = test_mine(mat_inp, mat_out, mat_out_total)
% test new

warning off
control=0;

a=[mat_inp];
b=[mat_out];
c=[mat_out_total];

input_binary=struct('filename',a,'control',control,'results_filename',b,'iter_num',2,'f',1.1,'max_num3',5,'th_last',0.4,'output_name',c);

% a is the cluster name
% b is the results file name (not all the variables)
% c is a filename to save all variables
% iter_num shows the circuit complexity

binary_regression4(input_binary)