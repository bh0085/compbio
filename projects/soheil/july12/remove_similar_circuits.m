function pf_new=remove_similar_circuits(pf,tf_sq)

% this function takes a profile pf and removes similar circuits withing
% threshold eps
eps=0.03;

if size(pf,1)==1 & length(pf{1,1})==0 | length([pf{:,1}])==0
    pf_new=pf;
% pf_new=[];
else
    row_size=size(pf,1);
    pf_new=cell(1,5);
    
    if length(pf{1,1})~=0
        pf_new(1,:)=pf(1,:);
    end
    if row_size>=2
        
        for i=2:row_size
            index_temp=size(pf_new,1)+1;
            sim_flag=0;
            y_i=read_module(pf{i,3},pf{i,4},tf_sq);
            for j=1:i-1
                
                if abs(pf{i,5}-pf{j,5})<=eps
                    
                    y_j=read_module(pf{j,3},pf{j,4},tf_sq);
                    if sum(sum(y_i~=y_j))<=eps.*length(y_j)
                        sim_flag=1;
                    end
                end
            end
            if sim_flag==0
                pf_new(index_temp,:)=pf(i,:);
            end
            
        end
    end
    
end


% if length([pf_new{:,3}])==0
%    pf_new=[]; 
% end