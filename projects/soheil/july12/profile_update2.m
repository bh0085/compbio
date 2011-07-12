function pf_new_th=profile_update2(pf,gene_sq,tf_sq,max_num,f)

% disp('inside profile update 2')
% pf{1,1}
% [pf{:,1}]

error_flag=0;
if length([pf{:,1}])==0
    
    cand_num=0;
    pf_new_th=cell(1,5);
else
    cand_num=size(pf,1);
    
    pf_new=pf;
    [ind,g_pattern]=gene_pattern(gene_sq);
    
    if length(ind)==0
%         disp('error in profile update, empty profile')
    end
    
    for i=1:cand_num
        i;
        for j=i:cand_num
            
            % and of these patterns
            [var_temp,code_temp]=my_and_modules(pf{i,3},pf{i,4},pf{j,3},pf{j,4});
            y=read_module(var_temp,code_temp,tf_sq);
            [q,p11,p22,p33,ind]=seq_qual2(y,gene_sq);
            
            if q>=f*max(pf{i,5},pf{j,5})
                
                if length(pf_new{1,1})~=0
                    index_temp=size(pf_new,1)+1; % row place to write
                else
                    index_temp=1;
                end
                pf_new{index_temp,1}=[ind];
                pf_new{index_temp,2}=[p11,p22,p33];
                pf_new{index_temp,3}=var_temp;
                pf_new{index_temp,4}=code_temp;
                pf_new{index_temp,5}=[q];
                
            end
            
            % or of these patterns
            [var_temp,code_temp]=my_or_modules(pf{i,3},pf{i,4},pf{j,3},pf{j,4});
            y=read_module(var_temp,code_temp,tf_sq);
            [q,p11,p22,p33,ind]=seq_qual2(y,gene_sq);
            
            if q>=f*max(pf{i,5},pf{j,5})
                
                
                if length(pf_new{1,1})~=0
                    index_temp=size(pf_new,1)+1; % row place to write
                else
                    index_temp=1;
                end
                pf_new{index_temp,1}=[ind];
                pf_new{index_temp,2}=[p11,p22,p33];
                pf_new{index_temp,3}=var_temp;
                pf_new{index_temp,4}=code_temp;
                pf_new{index_temp,5}=[q];
                
            end
%             % not(tf1)* tf2
%             %            disp('inside profile update')
%             %            size(pf{i,3})
%             %            size(pf{i,4})
%             [var_not,code_not]=my_not_modules(pf{i,3},pf{i,4});
%             [var_temp,code_temp]=my_and_modules(var_not,code_not,pf{j,3},pf{j,4});
%             y=read_module(var_temp,code_temp,tf_sq);
%             [q,p11,p22,p33,ind]=seq_qual2(y,gene_sq);
%             
%             if q>=f*max(pf{i,5},pf{j,5})
%                 
%                 
%                 if length(pf_new{1,1})~=0
%                     index_temp=size(pf_new,1)+1; % row place to write
%                 else
%                     index_temp=1;
%                 end
%                 pf_new{index_temp,1}=[ind];
%                 pf_new{index_temp,2}=[p11,p22,p33];
%                 pf_new{index_temp,3}=var_temp;
%                 pf_new{index_temp,4}=code_temp;
%                 pf_new{index_temp,5}=[q];
%                 
%             end
        end
    end
    
    pf_new_th=adjusted_profile_thresholding2(pf_new,max_num,tf_sq,gene_sq);
    
    %     if ind==0
    %
    %         pf0=profile_sorting(pf_new,[-1,0,1]);
    %         pf1=profile_sorting(pf_new,[1,0,-1]);
    %         pf2=profile_sorting(pf_new,[-1,0,0]);
    %         pf3=profile_sorting(pf_new,[1,0,0]);
    %         pf4=profile_sorting(pf_new,[0,0,1]);
    %         pf5=profile_sorting(pf_new,[0,0,-1]);
    %         pf6=profile_sorting(pf_new,[1,0,1]);
    %         pf7=profile_sorting(pf_new,[-1,0,-1]);
    %
    %
    %
    %     elseif ind==1 | ind==2
    %
    %         pf0=profile_sorting(pf_new,[-1,0,-2]);
    %         pf1=profile_sorting(pf_new,[1,0,-2]);
    %         pf2=profile_sorting(pf_new,[-1,1,-2]);
    %         pf3=profile_sorting(pf_new,[1,1,-2]);
    %         pf4=profile_sorting(pf_new,[0,1,-2]);
    %         pf5=profile_sorting(pf_new,[-1,-1,-2]);
    %         pf6=profile_sorting(pf_new,[0,-1,-2]);
    %         pf7=profile_sorting(pf_new,[1,-1,-2]);
    %
    %
    %
    %     elseif ind==3
    %         error_flag=1;
    %
    %     end
    %
    %     if error_flag==1
    %         pf_new_th=cell(1,5);
    %         disp('error!!! ind=3')
    %     else
    %
    %         [pf_new_0,th_used]=adjusted_profile_thresholding(pf0,max_num);
    %         [pf_new_1,th_used]=adjusted_profile_thresholding(pf1,max_num);
    %         [pf_new_2,th_used]=adjusted_profile_thresholding(pf2,max_num);
    %         [pf_new_3,th_used]=adjusted_profile_thresholding(pf3,max_num);
    %         [pf_new_4,th_used]=adjusted_profile_thresholding(pf4,max_num);
    %         [pf_new_5,th_used]=adjusted_profile_thresholding(pf5,max_num);
    %         [pf_new_6,th_used]=adjusted_profile_thresholding(pf6,max_num);
    %         [pf_new_7,th_used]=adjusted_profile_thresholding(pf7,max_num);
    %
    %
    %         %         disp('removing similarities')
    %         pf_new_0= remove_similar_circuits(pf_new_0,tf_sq);
    %         pf_new_1= remove_similar_circuits(pf_new_1,tf_sq);
    %         pf_new_2= remove_similar_circuits(pf_new_2,tf_sq);
    %         pf_new_3= remove_similar_circuits(pf_new_3,tf_sq);
    %         pf_new_4= remove_similar_circuits(pf_new_4,tf_sq);
    %         pf_new_5= remove_similar_circuits(pf_new_5,tf_sq);
    %         pf_new_6= remove_similar_circuits(pf_new_6,tf_sq);
    %         pf_new_7= remove_similar_circuits(pf_new_7,tf_sq);
    %
    %
    %
    %         pf_new_th= [pf_new_0;pf_new_1;pf_new_2;pf_new_3;pf_new_4;pf_new_5;pf_new_6;pf_new_7];
    %
    %
    %   disp('threshold inside the update')
    %         th_used
    
    
    
    %     end
    
end





