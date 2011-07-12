function [new_model,new_variables]=sing_check(model,variables)

beta=size(variables,1);
new_model=model;
new_variables=variables;


for i=1:beta
    if i<beta
        for j=i+1:beta
            
            if  sum(sum(variables{i,1}-variables{j,1}))==0
                new_model{1,j}=[];
                new_variables{j,1}=variables{j,1}.*0;
            end
        end
    end
end



