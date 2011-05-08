function [candidates_new,model_new,valid_genes]=update_candidates3(out_iter,candidates,error,error_en, gene,tf,net_prob,net_prob_sorted,net_sorted,gene_sim,f_mix,f_sim,model,k,check_names,trunc_value,valid_genes,beta)


[g_num,tf_num]=size(net_prob);

model_new=model;

corr_score=zeros(g_num,tf_num);
candidates_new=candidates;



if out_iter==1
    valid_genes=ones(1,g_num);
    for i=1:g_num
%         disp('inside update candid, iter=1')
        temp=find(net_prob_sorted(i,:)==0);
        if temp(1)>k
            candidates_new{i,1}=net_sorted(i,1:k);
        elseif temp(1)>beta
            candidates_new{i,1}=net_sorted(i,1:temp(1)-1);
        else
            valid_genes(i)=0;
%             disp('no TF for this gene')
        end
%        length(candidates_new{i,1}) 
        
    end
    
else
    
    for i=1:g_num
        if valid_genes(i)==1
            for j=1:tf_num
                temp=abs(corrcoef(tf{j,:},error{i,:}));
                corr_score(i,j)=temp(1,2);
            end
        end
    end
    
    
    new_net_prob=net_prob+f_mix.*corr_score;
    new_net_prob=new_net_prob.*check_names;
    %[new_net_sorted,new_net_prob_sorted]=net_sort(new_net_prob);
    
    for i=1:g_num
        if valid_genes(i)==1
            min_error_en=100000000;
            index=0;
            for j=1:g_num
                
                if gene_sim(i,j)==1 & error_en{j,1}(1,end)<min_error_en
                    min_error_en=error_en{j,1}(1,end);
                    index=j;
                end
            end
            
            if index~=0 & index~=i & error_en{i,1}(1,end)<min_error_en.*f_sim
                
                model_new(i,:)=model(index,:);
                temp=union(candidates{index,1},candidates{i,1});
                if length(temp)<=2*k
                    candidates_new{i,1}= temp;
                else
                    candidates_new{i,1}=[candidates{index,1}];
                    % we need to choose the best 2k-length(candidates{index,1})
                    % candidates from the old set
                    
                    temp_net_weight=zeros(1,k);
                    for j=1:k
                        temp_net_weight(j)=net_prob(i,candidates{i,1}(1,j));
                    end
                    [alaki,I]=sort(temp_net_weight,'descend');
                    
                    extra_cand=2*k-length(candidates{index,1});
                    if extra_cand>1
                        for j=1:extra_cand
                            candidates_new{i,1}=[candidates_new{i,1},candidates{i,1}(1,I(j))];
                        end
                    end
                    
                end
                
            else
%                             disp('MCMC')
                
                kk=length(candidates{i,1});
                temp=zeros(kk,tf_num);
                for j=1:kk
                    for l=1:tf_num
                        if ismember(l,candidates{i,1})==1
                            temp(j,l)=0;
                        else
%                             disp('inside update candid')
%                             i
%                             l
%                             kk
%                             candidates{i,1}(1,j)
                            temp(j,l)=max(new_net_prob(i,l)-new_net_prob(i,candidates{i,1}(1,j)),0);
                        end
                    end
                end
                
                a=sum(sum(temp));
                if a~=0
                    temp=temp./a; % MC prob matrix
                    
                    temp=mat_threshold(temp,trunc_value);
                    
                    [index1,index2]=rand_walk(temp);
                    
                    candidates_new{i,1}(1,index1)=index2;
                    
                end
                
                
            end
            
            
        end
    end
end