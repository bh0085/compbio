function pf_ss=profile_sorting_new(pf,pattern)

if length([pf{:,1}])~=0
    % this function gets a profile, separates the rows with the given pattern and then sorts them by their quality
    pf_s=cell(1,6);
    cand_num=size(pf,1);
    
    for i=1:cand_num
        if length(pf_s{1,1})~=0
            index_temp=size(pf_s,1)+1;
        else
            index_temp=1;
        end
        
        %    disp('inside profile sorting')
        %    size(pf{i,2})
        %    size(pattern)
%         disp('inside profile_sorting')
%         size(pf{i,2})
%         size(pattern)
        if length(pf{i,2})~=0
            
            if pf{i,2}==pattern
                pf_s(index_temp,:)=pf(i,:);
            end
        end
    end
    
    % now, we sort
    if length([pf_s{:,1}])~=0
        row_size=size(pf_s,1);
        quality_vector=zeros(row_size,1);
        
        for i=1:row_size
            quality_vector(i)=pf_s{i,5}*pf_s{i,6};
        end
        
        [alaki,I]=  sort(quality_vector,'descend');
        
        pf_ss=pf_s;
        for i=1:row_size
            pf_ss(i,:)=pf_s(I(i),:);
        end
        
    else
        pf_ss=pf_s;
    end
    
    
else
    pf_ss=cell(1,6);
end