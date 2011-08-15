% compute common modules accross different networks

clc
clear all
close all

a=load('rep_cell_denovo2');
a=a.rep_cell2;
size(a)

a=a(1:300,:);

% b=load('rep_cell_sush');
% sush=b.rep_cell2;
% size(sush)



% c=load('rep_cell_comb');
% comb=c.rep_cell2;
% size(comb)

b=load('rep_cell_inter_denovo2');
b=b.rep_cell2;
size(b)


ab=intersect_common_modules(a,b);
size(ab)

% denov_comb=intersect_common_modules(denovo,comb);
% size(denov_comb)
% 
% comb_sush=intersect_common_modules(comb,sush);
% size(comb_sush)
% 
% denov_sush_comb=intersect_common_modules(denov_sush,comb);
% size(denov_sush_comb)

length([ size(a,1),size(ab,1),size(b,1)])
vennX( [ size(a,1),size(ab,1),size(b,1)], .05 )

% saving data





% denov_sush=intersect_common_modules(denovo,sush);
% size(denov_sush)
% 
% denov_comb=intersect_common_modules(denovo,comb);
% size(denov_comb)
% 
% comb_sush=intersect_common_modules(comb,sush);
% size(comb_sush)
% 
% denov_sush_comb=intersect_common_modules(denov_sush,comb);
% size(denov_sush_comb)
% 
% length([ size(denovo,1),size(denov_sush,1),size(sush,1),size(comb_sush,1),size(comb,1),size(denov_comb,1),size(denov_sush_comb,1)])
% vennX( [ size(denovo,1),size(denov_sush,1),size(sush,1),size(comb_sush,1),size(comb,1),size(denov_comb,1),size(denov_sush_comb,1)], .05 )
% 
% % saving data
% 
% save('denovo','denovo');
% save('inter_cluster','sush');
% save('comb','comb');
% 
% save('denovo_inter','denov_sush');
% save('comb_inter','comb_sush');
% save('denovo_inter','denov_comb');
% 
% save('denovo_inter_comb','denov_sush_comb');


