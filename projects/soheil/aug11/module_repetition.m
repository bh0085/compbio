function result_cell=module_repetition(network_cell)

% result_cell is a cell with two colomns, first colomn is the module
% variables, the next colomn is the
% number of repetitions of the module

% circuit repetions should be added

row_size=size(network_cell,1);
if row_size==1 & length([network_cell{:,3}])==0
    result_cell=cell(1,2);
    
else
    result_cell=cell(1,2);
    for i=1:row_size
        done=0;
        for j=1:size(result_cell,1)
            if done==0 & length(result_cell{j,1})>=2 & length(network_cell{i,3})==length(result_cell{j,1})
                if sum(ismember(network_cell{i,3},result_cell{j,1}))==length(network_cell{i,3})
                    result_cell{j,2}=result_cell{j,2}+1;
                    done=1;
                end
            end
        end
        
        if done==0
            if length(result_cell{1,1})==0
                index_temp=1;
            else
                index_temp=size(result_cell,1)+1;
                
            end
            if length(network_cell{i,3})>=2
            result_cell{index_temp,1}=network_cell{i,3};
            result_cell{index_temp,2}=1;
            end
        end
    end
end

