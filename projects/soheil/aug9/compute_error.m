function [error,coefs]=compute_error(variables,gene)

gn=gene{1,1}';
beta=size(variables,1);
n=length(variables{1,1}(:));
variables_mat=zeros(beta,n);
for i=1:beta
    variables_mat(i,:)=variables{i,1};
end

stop=0;
for i=1:beta
    if variables_mat(i,:)==gene{1,1}
        stop=1;
    end
end

if stop==0
    
    % disp('inside compute error function')
    % size(variables_mat')
    size(gn);
    % variables_mat'
    % gn
    variables_mat=variables_mat';
    % disp('khaste!')
    % khasteh=gn-variables_mat(:,1)-3.*variables_mat(:,2)-variables_mat(:,3)-4.*variables_mat(:,4);
    % sum(sum(khasteh))
    
    new_variables_mat=[];
    indic=zeros(1,beta);
    for i=1:beta
        if variables_mat(:,i)~=0
            new_variables_mat=[new_variables_mat,variables_mat(:,i)];
            indic(i)=1;
            
        end
    end
    
    % disp('inside compute_error')
    % % size(variables_mat)
    % % % size(new_variables_mat)
    % indic;
    % variables_mat(1:5,:);
    % new_variables_mat
    % gn
    %
    %
    % new_variables_mat(1:5,:)
    % gn(1:5)
    [new_coefs,SIGMA,error] = mvregress(new_variables_mat,gn);
    
    coefs=zeros(1,beta);
    
    for i=1:beta
        if indic(i)~=0
            coefs(i)=new_coefs(sum(sum(indic(1:i))));
        end
    end
    
    
    error=error';
    ceofs=coefs';
    
else
    
    coefs=zeros(beta,1)+1000;
    error=zeros(1,n)+100000000;
    
end

