function [gene_sim,gene_sim_soft]=comp_gene_sim(net_prob,th_cor)

g_num=size(net_prob,1);
tf_num=size(net_prob,2);
gene_sim=zeros(g_num,g_num);


for i=1:g_num
    for j=1:g_num
        
%         a=net_prob(i,:)>th;%keeps the k_max values not k components
%         b=net_prob(j,:)>th;
% 
%             gene_sim(i,j)=sum(sum(a.*b+(1-a).*(1-b)))./tf_num;
            
        temp=abs(corrcoef(net_prob(i,:),net_prob(j,:)));
            gene_sim(i,j)=temp(1,2);    
 
        
    end
end

gene_sim_soft=gene_sim;

gene_sim=gene_sim>th_cor;