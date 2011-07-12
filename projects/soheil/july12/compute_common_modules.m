% compute common modules accross different networks

clc
clear all
close all

a=load('rep_cell_denovo');
denovo=a.rep_cell2;
size(denovo)

% b=load('rep_cell_sush');
% sush=b.rep_cell2;
% size(sush)



c=load('rep_cell_comb');
comb=c.rep_cell2;
size(comb)

b=load('rep_cell_inter');
sush=b.rep_cell2;
size(sush)

denov_sush=intersect_common_modules(denovo,sush);
size(denov_sush)

denov_comb=intersect_common_modules(denovo,comb);
size(denov_comb)

comb_sush=intersect_common_modules(comb,sush);
size(comb_sush)

denov_sush_comb=intersect_common_modules(denov_sush,comb);
size(denov_sush_comb)

length([ size(denovo,1),size(denov_sush,1),size(sush,1),size(comb_sush,1),size(comb,1),size(denov_comb,1),size(denov_sush_comb,1)])
vennX( [ size(denovo,1),size(denov_sush,1),size(sush,1),size(comb_sush,1),size(comb,1),size(denov_comb,1),size(denov_sush_comb,1)], .05 )

% saving data

save('denovo','denovo');
save('sush','sush');
save('comb','comb');

save('denov_sush','denov_sush');
save('comb_sush','comb_sush');
save('debov_comb','denov_comb');

save('denov_sush_comb','denov_sush_comb');


