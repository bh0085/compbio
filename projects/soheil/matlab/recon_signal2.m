function [sig,coefs]=recon_signal2(model,gene,tf)

% model is a 1*beta cell
%gene is a 1*1 cell
% tf is a tf_num*1 cell

n=length(gene{1,1});
beta=size(model,2);
variables=cell(beta,1);

% compute model variables

for i=1:beta
    if length(model{1,i})>=1
        temp=ones(1,n);
        for j=1:length(model{1,i})
            temp=temp.*tf{model{1,i}(1,j),1};
        end
        
        variables{i,1}=temp;
    else
        variables{i,1}=zeros(1,n);
    end
    
end

% for i=1:beta
%     
%     model{1,i}
% end
% 
% size(variables)
% for i=1:beta
%     
%     variables{i,1}
% end
% gene{1,1}

[error,coefs]=compute_error(variables,gene);

sig=zeros(1,n);
for i=1:beta
    sig=sig+coefs(i).*variables{i,1};
end


