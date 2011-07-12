function pf_th=adjusted_profile_thresholding2(pf,max_num,tf_sq,gene_sq)

% it threshold this profile to have at most num_max rows for each pattern

% disp('beginning of adjusted 2')
% 
% size(pf,1)
% length(pf{1,1})
% pf{1,1}


row_size=size(pf,1);
if row_size==1 & length(pf{1,1})==0 | length([pf{:,1}])==0
    pf_th=cell(1,5);
%  pf_th=[];
    
else
    
    [ind,g_pattern]=gene_pattern(gene_sq);
    error_flag=0;
    if ind==0
        
        pf0=profile_sorting(pf,[-1,0,1]);
        pf1=profile_sorting(pf,[1,0,-1]);
        pf2=profile_sorting(pf,[-1,0,0]);
        pf3=profile_sorting(pf,[1,0,0]);
        pf4=profile_sorting(pf,[0,0,1]);
        pf5=profile_sorting(pf,[0,0,-1]);
        pf6=profile_sorting(pf,[1,0,1]);
        pf7=profile_sorting(pf,[-1,0,-1]);
        
        
        
    elseif ind==1 | ind==2
        
        pf0=profile_sorting(pf,[-1,0,-2]);
        pf1=profile_sorting(pf,[1,0,-2]);
        pf2=profile_sorting(pf,[-1,1,-2]);
        pf3=profile_sorting(pf,[1,1,-2]);
        pf4=profile_sorting(pf,[0,1,-2]);
        pf5=profile_sorting(pf,[-1,-1,-2]);
        pf6=profile_sorting(pf,[0,-1,-2]);
        pf7=profile_sorting(pf,[1,-1,-2]);
        
        
        
    elseif ind==3
        error_flag=1;
        
    end
    
    % ploting the heatmap
%     mat=[gene_sq
%         read_module(pf0{1,3},pf0{1,4},tf_sq)
%         read_module(pf1{1,3},pf1{1,4},tf_sq)
%         read_module(pf2{1,3},pf2{1,4},tf_sq)
%         read_module(pf3{1,3},pf3{1,4},tf_sq)
%         read_module(pf4{1,3},pf4{1,4},tf_sq)
%         read_module(pf5{1,3},pf5{1,4},tf_sq)
%         read_module(pf6{1,3},pf6{1,4},tf_sq)
%         read_module(pf7{1,3},pf7{1,4},tf_sq)];
%     
% heatmap(mat,'RowLabels',1:size(mat,1),'ColumnLabels',1:size(mat,2),'Colormap',redbluecmap(11))

    
    if error_flag==1
        pf_th=cell(1,5);
%         disp('error!!! ind=3')
        
    else
%         disp('inside adjusted 2')
%         ind
%         size(pf0)
%         pf0
        [pf_new_0,th_used]=adjusted_profile_thresholding(pf0,max_num);
        [pf_new_1,th_used]=adjusted_profile_thresholding(pf1,max_num);
        [pf_new_2,th_used]=adjusted_profile_thresholding(pf2,max_num);
        [pf_new_3,th_used]=adjusted_profile_thresholding(pf3,max_num);
        [pf_new_4,th_used]=adjusted_profile_thresholding(pf4,max_num);
        [pf_new_5,th_used]=adjusted_profile_thresholding(pf5,max_num);
        [pf_new_6,th_used]=adjusted_profile_thresholding(pf6,max_num);
        [pf_new_7,th_used]=adjusted_profile_thresholding(pf7,max_num);
        
        
        %         disp('removing similarities')
        pf_new_0= remove_similar_circuits(pf_new_0,tf_sq);
        pf_new_1= remove_similar_circuits(pf_new_1,tf_sq);
        pf_new_2= remove_similar_circuits(pf_new_2,tf_sq);
        pf_new_3= remove_similar_circuits(pf_new_3,tf_sq);
        pf_new_4= remove_similar_circuits(pf_new_4,tf_sq);
        pf_new_5= remove_similar_circuits(pf_new_5,tf_sq);
        pf_new_6= remove_similar_circuits(pf_new_6,tf_sq);
        pf_new_7= remove_similar_circuits(pf_new_7,tf_sq);
        

        
        
        
        pf_th= [pf_new_0;pf_new_1;pf_new_2;pf_new_3;pf_new_4;pf_new_5;pf_new_6;pf_new_7];
        

        
        
        
    end
    
    
end

