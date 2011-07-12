function result_cell=intersect_common_modules(rep_cell1,rep_cell2)

% rep_cell1 and rep_cell2 are cells of some rows and 2 colomns
%the first colomn is the module, the second colomn is the number of
%repetions of that module
% result is a cell with some rows and two colomns, the first colomn is a
% common module and the second colomn is a vector of repetion numbers

n1=size(rep_cell1,1);
n2=size(rep_cell2,1);
result_cell=cell(1,2);

if n1>=1 & n2>=1 
    for i=1:n1
        done=0;
        for j=1:n2
            if done==0 & length(rep_cell1{i,1})==length(rep_cell2{j,1}) & sum(ismember(rep_cell1{i,1},rep_cell2{j,1}))==length(rep_cell2{j,1})
                if length(result_cell{1,1})~=0
                    index_temp=size(result_cell,1)+1;
                else
                    index_temp=1;
                end
                
                result_cell{index_temp,1}=rep_cell1{i,1};
                result_cell{index_temp,2}=[rep_cell1{i,2},rep_cell2{j,2}];
                done=1;
            end
        end
    end
    
else
    
    disp('empty inputs in intersect_common_modules!')
end
