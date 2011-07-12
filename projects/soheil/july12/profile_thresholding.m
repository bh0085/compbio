function pf_th=profile_thresholding(pf,th)

pf_th=cell(1,5);
row_size=size(pf,1);
if length([pf{:,1}])==0
    row_size=0;
end

if row_size~=0
    for i=1:row_size
        if length(pf_th{1,1})~=0
            index_temp=size(pf_th,1)+1;
        else
            index_temp=1;
        end
        
        if pf{i,5}>=th
            pf_th(index_temp,:)=pf(i,:);
        end
        
    end
    
    
end

% if length([pf_th{:,3}])==0
%    pf_th=[]; 
% end