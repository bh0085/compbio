function [network_cell,infered_network,quality_test]=binary_regression5(input)

filename=input.filename;

results_filename=input.results_filename;

output_filename=input.output_name;

control=input.control;

prior=input.prior;

% network_name=input.network_name;

% it adjusts thresholding level to have around this many candidates
% max_num1=input.max_num1;
% max_num2=input.max_num2;
max_num3=input.max_num3;
% max_num_last=input.max_num_last;
th_last=input.th_last;


iter_num=input.iter_num;
f=input.f;

% [gene_tot,tf_tot,net_prob,gene_names,tf_names]=prepare_data2(filename);
[gene_tot,tf_tot,net_prob,gene_names,tf_names]=prepare_data6(filename,control);


g_num=size(gene_tot,1);
tf_num=size(tf_tot,1);

check_names=ones(g_num,tf_num);
% zero in this matrix means they are the same genes
for i=1:g_num
    for j=1:tf_num
        temp_names=gene_names{i,:}==tf_names{j,:};
        if ismember(0,temp_names)==0
            check_names(i,j)=0;
        end
    end
end


% gene_test,tf_test

len=length(gene_tot{1,1});

test_times=1:4:len;
train_times=[];
for i=1:len
    if ismember(i,test_times)~=1
        train_times=[train_times,i];
    end
end

% test_times=floor(len*rand(1,len_test));
gene_test=cell(g_num,1);
tf_test=cell(tf_num,1);

for i=1:g_num
    gene_test{i,1}=gene_tot{i,1}(test_times);
end

for i=1:tf_num
    tf_test{i,1}=tf_tot{i,1}(test_times);
end

%********************

gene=cell(g_num,1);
tf=cell(tf_num,1);

for i=1:g_num
    gene{i,1}=gene_tot{i,1}(train_times);
end

for i=1:tf_num
    tf{i,1}=tf_tot{i,1}(train_times);
end

%******************
% results
network_cell=cell(1,6);
% gene_name|not_flag=0 (pattern)
% not_flag=1(not of the pattern)|circuit variables|circuit module|quality
% 6th cell is the prior of that module
%*****************


valid_genes=zeros(1,g_num);


for i=1:g_num
    if sum(net_prob(i,:)>0)>=2
        valid_genes(i)=1;
    end
end

% num_valid_genes=sum(valid_genes);


% [gene_s4,I]=sort(gene{4,1});
%     gene_sq4=quantize_gaussian2(gene_s4,1);
% 
% 
% [gene_s10,I]=sort(gene{10,1});
%     gene_sq10=quantize_gaussian2(gene_s10,1);
%     
%     gene_zero=gene_sq10*0;
%     
%     mat=[gene_sq10;gene_zero;gene_zero;gene_zero];
% heatmap(mat,'Colormap',redbluecmap(11))
% 
%  mat2=[gene_sq4;gene_zero;gene_zero;gene_zero];
% heatmap(mat2,'Colormap',redbluecmap(11))

for i=1:g_num
    if valid_genes(i)==1
    i
    valid_tf=zeros(1,tf_num);
    for j=1:tf_num
        if net_prob(i,j)>0 & check_names(i,j)==1
        valid_tf(j)=1;
        end
    end
    
    % sort and quantize
    [gene_s,I]=sort(gene{i,1});
    gene_sq=quantize_gaussian2(gene_s,1);
    
%     heatmap(gene_sq,'Colormap',redbluecmap(11))

    
    tf_sq=zeros(tf_num,length(gene_s));
    
    for j=1:tf_num
        if valid_tf(j)==1
        tf_sq(j,:)=quantize_gaussian2(tf{j,1}(I),1);
        end
    end
    
    
    
    % determining the type of gene pattern
    
    [ind,g_pattern]=gene_pattern(gene_sq);
    if g_pattern==[-2,-2,-2]
    disp('constant gene profile')
    end
    % ind is the kind of pattern that we have
    
    
    pf=cell(1,5); % first cell is the ind
    %second cell is the profile pattern
    % third var names
    % fourth: the binary function
    % fifth: the quality
    
    % initialization with ind. TF's
    
    for j=1:tf_num
        if valid_tf(j)==1
            [q,p1,p2,p3,ind]=seq_qual2(tf_sq(j,:),gene_sq);
            
            
            if length(pf{1,1})~=0
                index_temp=size(pf,1)+1; % row place to write
            else
                index_temp=1;
            end
            if q>0
            
            pf{index_temp,1}=[ind];
            pf{index_temp,2}=[p1,p2,p3];
            pf{index_temp,3}=[j];
            pf{index_temp,4}=[j];
            pf{index_temp,5}=[q];
            pf{index_temp,6}=compute_prior(prior, pf{index_temp,4});
            end
            
        end
    end
    
%     disp('profile size after ind')
%     g_pattern
% 
%     size(pf,1)
%     [pf{:,5}]'
%     [pf{:,2}]'
%     plot_profile(gene_sq,tf_sq,pf,5);

    pf=adjusted_profile_thresholding2_new(pf,max_num3,tf_sq,gene_sq);
    var_sing=profile_variables(pf);
    
%     disp('profile size after ind adjust')
%     size(pf,1)
%     [pf{:,5}]'
%     [pf{:,2}]'
%     ploting initial TF's
    
%         disp('after initialization')
%         pf_s1=profile_sorting(pf,g_pattern);
%         pf_s2=profile_sorting(pf,my_not(g_pattern));
%     
%         pf_init=[pf_s1;pf_s2];
%         %     size(pf_s1_th,1)==size(pf_s1,1)
%     
%     
%         [pf_init{:,5}]'
%     
%     
%         plot_profile(gene_sq,tf_sq,pf_init,5);
    %
    %
    
    
    %         heatmap_test1=[gene_sq;tf_sq(17,:)];
    %
    %             heatmap(heatmap_test1,'RowLabels',1:size(heatmap_test1,1),'ColumnLabels',1:size(heatmap_test1,2),'Colormap',redbluecmap(11))
    
    
    
    
    
    % initialization with TF modules without the high quality TF's
    
    for j1=1:tf_num
        for j2=j1+1:tf_num
            if valid_tf(j1)==1 & valid_tf(j2)==1 & ismember(j1,var_sing)==0 & ismember(j2,var_sing)==0
                
                
                
                % and
                y=my_and(tf_sq(j1,:),tf_sq(j2,:));
                [q,p11,p22,p33,ind]=seq_qual2(y,gene_sq);
                
                
                
                if length(pf{1,1})~=0
                    
                    index_temp=size(pf,1)+1; % row place to write
                else
                    index_temp=1;
                end
                if q>0
                pf{index_temp,1}=[ind];
                pf{index_temp,2}=[p11,p22,p33];
                pf{index_temp,3}=[j1,j2];
                pf{index_temp,4}=[j1,j2];
                pf{index_temp,5}=[q];
                pf{index_temp,6}=compute_prior(prior, pf{index_temp,4});
                end
                % or
                y=my_or(tf_sq(j1,:),tf_sq(j2,:));
                [q,p11,p22,p33,ind]=seq_qual2(y,gene_sq);
                
                
                if length(pf{1,1})~=0
                    index_temp=size(pf,1)+1; % row place to write
                else
                    index_temp=1;
                end
                if q>0
                pf{index_temp,1}=[ind];
                pf{index_temp,2}=[p11,p22,p33];
                pf{index_temp,3}=[j1,j2];
                pf{index_temp,4}=[j1,-1,j2];
                pf{index_temp,5}=[q];
                pf{index_temp,6}=compute_prior(prior, pf{index_temp,4});
                end
                
                
            end
        end
    end
    
%     disp('profile size after modules')
%     size(pf,1)
%     
%     
    pf=adjusted_profile_thresholding2_new(pf,max_num3,tf_sq,gene_sq);
    %     
%     disp('profile size after modules adjust')
%         plot_profile(gene_sq,tf_sq,pf,10);
% 
%     size(pf,1)
%     [pf{:,5}]
%     [pf{:,2}]
    % infering the function
    if iter_num~=0
        for iter=1:iter_num
    iter;        
            pf=profile_update3(pf,gene_sq,tf_sq,max_num3,f,prior);
            
%    disp('profile size after update')
%     size(pf,1)
%     [pf{:,5}]  
%     [pf{:,2}]
%     plot_profile(gene_sq,tf_sq,pf,10);

    
        end
    end
    pf=remove_empty_rows_new(pf);
    % 
%     disp('last profile')
%     pf{:,:}
    % pick the highest quality pf0
    
    pf_num=size(pf,1);
    
    if pf_num==0
%         disp('no infered profile')
    else
        pf_last=profile_sorting_new(pf,g_pattern);
        pf_last_not=profile_sorting_new(pf,my_not(g_pattern));
    end
    
    pf_last_2=pf_last_not;
    %     if size(pf_last_not,1)>=1 & length(pf_last_not{1,1})~=0
    %         for j=1:size(pf_last_not,1)
    %             pf_last_2{j,2}=my_not(pf_last_not{j,2});
    %             [pf_last_2{j,3},pf_last_2{j,4}]=my_not_modules(pf_last_not{j,3},pf_last_not{j,4});
    %         end
    %
    %     end
    
    pf_circuit=[pf_last;pf_last_2];
    
    pf_circuit= remove_similar_circuits_new(pf_circuit,tf_sq);
    pf_circuit2=profile_thresholding_new(pf_circuit,th_last);
    
    
%   if no model, choose the best existent one  
    if length([pf_circuit2{:,3}])~=0
    pf_circuit=pf_circuit2;
    
    elseif length([pf_circuit{:,3}])==0
        
        pf_circuit=profile_thresholding_new(pf,th_last);        
    elseif length([pf_circuit{:,3}])==0
        
        pf_circuit=profile_thresholding_new(pf,th_last/2);
    end
    
    
    
    
    
%     [pf_circuit,th_used2]=adjusted_profile_thresholding(pf_circuit,max_num_last);
    
%         disp('final circuit size')
%     size(pf_circuit,1)
%     [pf_circuit{:,:}]
    
    %     alaki1=[pf_circuit{:,5}]'
    %
%         plot_profile(gene_sq,tf_sq,pf_circuit,10);
    %
    %****************
    % saving the results for gene i
    
    
    
    if size(pf_circuit,1)>=1
        for j=1:size(pf_circuit,1);
            if length(pf_circuit{j,3})~=0
                
                row_size=size(network_cell,1);
                if row_size==1 & length(network_cell{1,1})==0
                    index_temp=1;
                else
                    index_temp=row_size+1;
                end
                
                if pf_circuit{j,2}==g_pattern
                    not_flag=0;
                elseif pf_circuit{j,2}==my_not(g_pattern)
                    not_flag=1;
                else
                    not_flag=2; % not matching the pattern
                end
                
                network_cell{index_temp,1}=i; % gene name
                network_cell{index_temp,2}=not_flag;
                network_cell{index_temp,3}=pf_circuit{j,3};
                network_cell{index_temp,4}=pf_circuit{j,4};
                network_cell{index_temp,5}=pf_circuit{j,5};
                network_cell{index_temp,6}=pf_circuit{j,6};
    
            end
        end
    end
    
end   
end
%********************************
infered_network=zeros(g_num,tf_num);
row_size=size(network_cell,1);
if row_size>=1
    for i=1:row_size
        gene_id=network_cell{i,1};
        var=network_cell{i,3};
        if length(var)>=1
            for j=1:length(var)
                infered_network(gene_id,var(j))=1;
            end
        end
    end
    
end

%****************************************
% performance on the test set
gene_test_q=zeros(g_num,length(gene_test{1,1}));

for i=1:g_num
    gene_test_q(i,:)=quantize_gaussian2(gene_test{i,1},1);
end

tf_test_q=zeros(tf_num,length(tf_test{1,1}));

for i=1:tf_num
    tf_test_q(i,:)=quantize_gaussian2(tf_test{i,1},1);
end

quality_test=zeros(g_num,1);
row_size=size(network_cell,1);

if row_size>=1
    for i=1:row_size
        gene_id=network_cell{i,1};
        var=network_cell{i,3};
        cir=network_cell{i,4};
        if length(var)~=0
            y=read_module(var,cir,tf_test_q);
            quality_temp=sum(sum(y==gene_test_q(gene_id,:)))/length(gene_test_q(gene_id,:));
            quality_test(gene_id,1)=max(quality_test(gene_id,1),quality_temp);
        end
    end
    
end

% % error correcting code?
% 
% 
% quality_test2=zeros(g_num,1);
% quality_test3=zeros(g_num,1);
% m=length(tf_test{1,1});
% 
% for i=1:g_num
%     cir_index=[];
%     if row_size>=1
%         for j=1:row_size
%             if network_cell{j,1}==i
%                 cir_index=[cir_index,j];
%             end
%         end
%     end
%     
%     
%     len=length(cir_index);
%     if len>=1
%         y=zeros(len,m+1); % last colomn is the quality of each row
%         for j=1:len
%             y(j,1:end-1)=read_module(network_cell{cir_index(j),3},network_cell{cir_index(j),4},tf_test_q);
%             y(j,end)=network_cell{cir_index(j),5};
%         end
%         
% 
% index_temp=find(y(:,end)==max(y(:,end)));        
% index_temp=index_temp(1);
% if length(index_temp)~=0
% quality_test3(i)=sum(sum(y(index_temp,1:end-1)==gene_test_q(i,:)))/m;
% else
%     quality_test3(i)=0;
% end
% 
%         
%         y_final=zeros(1,m);
%         for j=1:m
%             a=sum(sum(y(:,j)==0));
%             b=sum(sum(y(:,j)==1));
%             c=sum(sum(y(:,j)==-1));
%             
%             index_a=find(y(:,j)==0);
%             f_a=length(index_a);
%             q_a=y(index_a,m);
%             if f_a~=0
%             q_a=sum(sum(q_a))/f_a;
%             else
%              q_a=0;   
%             end
%             
%             
%             index_b=find(y(:,j)==1);
%             f_b=length(index_b);
%             q_b=y(index_b,m);
%             if f_b~=0
%             q_b=sum(sum(q_b))/f_b;
%             else
%              q_b=0;   
%             end
%             
%             index_c=find(y(:,j)==-1);
%             f_c=length(index_c);
%             q_c=y(index_c,m);
%             if f_c~=0
%             q_c=sum(sum(q_c))/f_c;
%             else
%              q_c=0;   
%             end
%             
% %             if (f_a*q_a)> (f_b*q_b) & (f_a*q_a)> (f_c*q_c)
% %                 y_final(j)=0;
% %             elseif (f_b*q_b)> (f_a*q_a) & (f_b*q_b)> (f_c*q_c)
% %                 y_final(j)=1;
% %             elseif (f_c*q_c)> (f_b*q_b) & (f_c*q_c)> (f_a*q_a)
% %                 y_final(j)=-1;
% %                 
% %             else
% %                 index_temp=find(y(:,end)==max(y(:,end)));
% %                 index_temp=index_temp(1);
% %                 y_final(j)=y(index_temp,j); 
% %             end
% 
%             [x,f]=mode(y(:,j));
%             
% 
%             if (f-a)==0 & (f-b)(f-c)~=0;
%                 y_final(j)=0;
%             elseif (f-b)==0 & (f-a)(f-c)~=0;
%                 y_final(j)=1;
%             elseif (f-c)==0 & (f-b)(f-a)~=0;
%                 y_final(j)=-1;
%             else
%                 
%                 index_temp=find(y(:,end)==max(y(:,end)));
%                 index_temp=index_temp(1);
%                 y_final(j)=y(index_temp,j);
%             end
% 
%         end
%         
%         
%         quality_test2(i)=sum(sum(y_final==gene_test_q(i,:)))/m;
% 
%         
%     else
%         y_final=zeros(1,m)-2;
%         quality_test2(i)=0;
%     end
%     
% end
% 
% 
% 
% 
% mean(quality_test)
% mean(quality_test2)
% mean(quality_test3)
% 
% 
% 
% figure
% plot(quality_test)



% figure
% hist(quality_test)
% 
% 
% % figure
% % plot(quality_test2)
% figure
% hist(quality_test2)
% 
% % figure
% % plot(quality_test3)
% figure
% hist(quality_test3)
% 

%     %*********************************************
%     % test this model over unused data (test data)
%     
%     
%     [gene_test_s,I]=sort(gene_test{i,1});
%     gene_test_sq=quantize_gaussian2(gene_test_s,1);
%     
%     tf_test_sq=zeros(tf_num,length(gene_test_s));
%     
%     for j=1:tf_num
%         tf_test_sq(j,:)=quantize_gaussian2(tf_test{j,1}(I),1);
%     end
%     
%     
%     %    [ind_test,g_test_pattern]=gene_pattern(gene_test_sq);
%     
%     
%     cir_num=size(pf_circuit,1);
%     if cir_num==1 & length(pf_circuit{1,1})==0
%         disp('no circuit to evaluate')
%     else
%         quality_vector=zeros(cir_num,1);
%         index_test=[];
%         for j=1:cir_num
%             %             disp('inside binary_regression')
%             %             size(pf_circuit{j,2})
%             %             size(g_pattern)
%             if length(pf_circuit{j,2})~=0
%                 if pf_circuit{j,2}==g_pattern
%                     y=read_module(pf_circuit{j,3},pf_circuit{j,4},tf_test_sq);
%                     quality_vector(j)=sum(sum(y==gene_test_sq))/length(gene_test_sq);
%                 elseif pf_circuit{j,2}==my_not(g_pattern)
%                     y=my_not(read_module(pf_circuit{j,3},pf_circuit{j,4},tf_test_sq));
%                     quality_vector(j)=sum(sum(y==gene_test_sq))/length(gene_test_sq);
%                     
%                 end
%             else
%                 quality_vector(j)=0;
%             end
%         end
%         
%     end
%     quality_vector
%     
%     % displaying test results
%     
%     plot_profile(gene_test_sq,tf_test_sq,pf_circuit,10);
%     max(quality_vector)
%     quality_test(i)= max(quality_vector);
%     
% 
% 
% 
% 
% quality_test

cluster_size=length(gene_tot{1,1});

save(output_filename)
save(results_filename,'network_cell','infered_network','quality_test','cluster_size')












