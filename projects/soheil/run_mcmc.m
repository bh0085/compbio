function run_mcmc(mat_inp, mat_out)
%Takes filenames of:
% input: a file that must currently exist on the filesystem
%        containing a single matlab struct "parameters"
%
% output:a file that need not yet exist to which we will write
%        execution results for mat_inp
%
input = open(mat_inp)

%output =....



%% test 7
%clc
%clear all
%close all
%warning off
%
%input=struct('filename','expression_c4d_n4_tt_4.mat','out_iter_num',2,'in_iter_num',5,'k',6,'beta',4,'f_mix',1,'f_sim',0.9,'f_en_in',0.95,'f_en_out',0.95,'th_cor',0.5,'trunc_value',3,'degree_bound',3,'output_name','result');
%
%%out_iter_num=50
%%in_iter_num=1000
%
%%k=4:10
%% beta=4:k beta is less or equal to k
%
%%f_mix=0.5:0.1:2
%%f_sim= 0.5:0.05:0.95
%%f_en_in=0.5:0.05:1
%%f_en_out=0.5:0.05:1
%%th_cor=0.4:0.05:0.8
%%trunc_value=2:5
%%degree_bound=3:6
%
%MCMC_noonlinear(input);


filename=input.filename;
out_iter_num=input.out_iter_num;
in_iter_num=input.in_iter_num;
k=input.k; % 2k is the bound on the number of used candidates 

beta=input.beta; % model order always less or equal to k
f_mix=input.f_mix;% scale for the candidate updating by using correlation coefs alpha zero means no error correlation effect in candidate updating
f_sim=input.f_sim; % gene similarity comparsion factor (smaller this ratio, harder to use similarity)
f_en_in=input.f_en_in; %improvement factor
f_en_out=input.f_en_out;

th_cor=input.th_cor;
trunc_value=input.trunc_value; % inside the update_candidate, in MCMC part, it looks at the truc_value edges with the highest probabilities
degree_bound=input.degree_bound;

b=beta;


[gene_tot,tf_tot,net_prob,gene_names,tf_names]=prepare_data2(filename);

g_num=size(gene_tot,1);
tf_num=size(tf_tot,1);

check_names=ones(g_num,tf_num);
% one in this matrix means they are the same genes
for i=1:g_num
    for j=1:tf_num
        temp_names=gene_names{i,:}==tf_names{j,:};
       if ismember(0,temp_names)==0
          check_names(i,j)=0;
       end
    end
end

net_prob=net_prob.*check_names;

% zero mean
for i=1:g_num
    gene_tot{i,1}=gene_tot{i,1}-mean(gene_tot{i,1});
end

for i=1:tf_num
    tf_tot{i,1}=tf_tot{i,1}-mean(tf_tot{i,1});
end

% unite variance

for i=1:g_num
    gene_tot{i,1}=gene_tot{i,1}./sqrt(mean(gene_tot{i,1}.^2));
end

for i=1:tf_num
    tf_tot{i,1}=tf_tot{i,1}./sqrt(mean(tf_tot{i,1}.^2));
end


% gene_test,tf_test

len=length(gene_tot{1,1});
len_test=floor(len/4);
test_times=1:floor(len/4):len;
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

gene=cell(g_num,1);
tf=cell(tf_num,1);

for i=1:g_num
    gene{i,1}=gene_tot{i,1}(train_times);
end

for i=1:tf_num
    tf{i,1}=tf_tot{i,1}(train_times);
end


[net_sorted,net_prob_sorted]=net_sort(net_prob);

gene_sim=zeros(g_num,g_num); % 0-1 matrix
[gene_sim,gene_sim_soft]=comp_gene_sim(net_prob,th_cor);

candidates=cell(g_num,1);
candidates_new=cell(g_num,1);

error_new=cell(g_num,1); % error term of gene expression before the update
error=cell(g_num,1);

error_en=cell(g_num,1);
error_en_new=cell(g_num,1);

model_storage=cell(g_num,out_iter_num*in_iter_num);
candid_storage=cell(g_num,out_iter_num);


for i=1:g_num
    error_en{i,1}=100000000;
    error_en_new{i,1}=100000000;
end


model=cell(g_num,beta);
model_linear=cell(g_num,beta);
model_new=cell(g_num,beta);

valid_genes=[];

for out_iter=1:out_iter_num
   [candidates_new,model,valid_genes]= update_candidates3(out_iter,candidates,error,error_en, gene,tf,net_prob,net_prob_sorted,net_sorted,gene_sim,f_mix,f_sim,model,k,check_names,trunc_value,valid_genes,beta);

   
   for i=1:g_num
       candid_storage{i,out_iter}=candidates_new(i,1);
   end
   
   
   
    for i=1:g_num
        % choose the model
        i;
        if valid_genes(i)==1

        length(candidates_new{i,1});
        variables=cell(beta,1);
        
        for in_iter=1:in_iter_num
            in_iter;
            stop=0; % when we don't have any direction to go in inside loop it becomes 1
            if stop==0
                

                [new_model,new_variables,stop,model_linear_i]=update_model3(in_iter,out_iter, model(i,:),error_new(i,1),candidates_new{i,1},tf,gene(i,1),variables,b,f_en_in,degree_bound);
                
                if in_iter==1 & out_iter==1
                    model_linear(i,:)=model_linear_i;
                end

                [new_model,new_variables]=sing_check(new_model,new_variables); % check the singularity of the variables
                [prid_error,coefs]=compute_error(new_variables,gene(i,1));
                variables=new_variables;
                model(i,:)=new_model;
                error_new{i,1}=prid_error;
                error_en_new{i,1}=[error_en_new{i,1},sum(sum(prid_error.^2))];
                model_storage{i,(out_iter-1)*in_iter_num+in_iter}=new_model;
                
            end
            
            
        end
    
    
    
    % to accept or reject these candidates
    % compute the error energy for each gene
    %if less than previos one accept new candidates and replace error term
    % otherwise stich with previous candidates
    
    
  
        if error_en_new{i,1}(1,end)<f_en_out.*error_en{i,1}(1,end)
            candidates{i,1}=candidates_new{i,1};
            error_en{i,1}=[error_en{i,1},error_en_new{i,1}(1,end)];
            error(i,1)=error_new(i,1);
        end
    end
    end 
    
end






% display results for gene s


len=length(error_en_new{1,1})-1;
total_error=zeros(g_num,len);

for i=1:g_num
    if valid_genes(i)==1 
    total_error(i,:)=error_en_new{i,1}(2:end);
    end
end


%heatmap(total_error,'RowLabels',1:size(total_error,1),'ColumnLabels',1:size(total_error,2),'Colormap',redbluecmap(11))

ratio=zeros(1,g_num);
for i=1:g_num
    if valid_genes(i)==1
    ratio(i)=total_error(i,end)/total_error(i,1);
    end
end


error_alaki=zeros(g_num,1);



% validation of our model for the test data
error_test=zeros(g_num,2);
coefs_dic_linear=cell(g_num,1);
coefs_dic_nonlinear=cell(g_num,1);
for i=1:g_num
    if valid_genes(i)==1
        
        [alaki,coefs_linear]=recon_signal2(model_linear(i,:),gene(i,:),tf);
        [alaki,coefs_nonlinear]=recon_signal2(model(i,:),gene(i,:),tf);
        
        coefs_dic_linear{i,1}=coefs_linear;
        coefs_dic_nonlinear{i,1}=coefs_nonlinear;
        
        sig_rec_linear=recon(model_linear(i,:),gene_test(i,:),tf_test,coefs_linear);
        sig_rec_nonlinear=recon(model(i,:),gene_test(i,:),tf_test,coefs_nonlinear);
        
        
        error_linear=gene_test{i,:}-sig_rec_linear;
        error_nonlinear=gene_test{i,:}-sig_rec_nonlinear;
        error_test(i,1)=sum(sum(error_linear.^2));
        error_test(i,2)=sum(sum(error_nonlinear.^2));
    end
end


ratio_valid=zeros(1,g_num);
for i=1:g_num
    if error_test(i,1)~=0 & valid_genes(i)==1
        ratio_valid(i)=error_test(i,2)./error_test(i,1);
    end
end


b1=sum(sum(ratio_valid==0));
b2=sum(sum(ratio_valid==1));
b3=sum(sum(ratio_valid>1));

total_num=g_num-b1;
improv_ratio=(g_num-b1-b2-b3)/total_num;
stay_same=b2/total_num;
get_worst_ratio=1-improv_ratio-stay_same;


save(mat_output)









