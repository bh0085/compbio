function [new_model,new_variables,stop,model_linear]=update_model3(in_iter,out_iter, model,error,candid,tf,gene,variables,b,c,degree_bound)


beta=size(model,2);
score=zeros(beta,1);
kk=length(candid); % note: this kk is different than the global defined k
trans_mat=zeros(beta,2*kk);
ss=zeros(beta,1);
% b the number of considered model variables, b \leq beta
%variables=cell(beta,1);
% stop related to linear model saving
new_variables=cell(beta,1);

model_linear=cell(1,beta);


if in_iter==1 & out_iter==1
    
    for i=1:beta
        if i<=length(candid)
        new_model{1,i}=[candid(i)];
        variables(i,1)=tf(candid(i),1);
        new_variables(i,1)=variables(i,1);
        else
          variables{i,1}=tf{1,1}.*0;
          new_variables(i,1)=variables(i,1);
        end
    end
    
    model_linear=new_model;
    
    stop=0;
    
elseif in_iter==1
    
    % update new model and new variables by output of the outer loop
    % here you need to compute variables again
    new_model=model;
    
    for i=1:beta
        
        if length(model{1,i})~=0
            
            temp_var=ones(1,length(tf{1,1}));
            for j=1:length(model{1,i})
                temp_var=temp_var.*tf{model{1,i}(1,j),1};
            end
        else
            temp_var=zeros(1,length(tf{1,1}));
        end
        new_variables{i,1}=temp_var;
    end
    % note:
    % here we could compare the result of this with the linear one to see if it is acceptable or not
    
    stop=0;
else
    
    error_energy=zeros(beta,2*kk)+100000000;
    
    for i=1:beta
        temp1=(corrcoef(error{1,1},variables{i,1}));
        temp2=temp1(1,2);
        score(i,1)=abs(temp2);
    end
    
    [sorted_score,I]=sort(score,'descend');
    
    
    for i=1:b
        for j=1:kk
            % compute the error if we change Z(I(i)) to Z(I(i)).*candid(j)
            if length(model{1,I(i)})<degree_bound
            new_model=model;
            new_variables=variables;
            new_model{1,I(i)}=[model{1,I(i)},candid(j)];
            new_variables{I(i),1}=variables{I(i),1}.*tf{candid(j),1};
            
%             disp('model inside the model update')
%             
%             
%             for kkk=1:beta
%                 new_model{1,kkk}
%             end
            [new_model,new_variables]=sing_check(new_model,new_variables);
%             disp('model inside the model update11')
%             
%             
%             for kkk=1:beta
%                 new_model{1,kkk}
%             end
            [prid_error,coefs]=compute_error(new_variables,gene);
            error_energy(I(i),j)=sum(sum(prid_error.^2));
            end
        end
        
        for j=1:kk
            
            new_model=model;
            new_variables=variables;
            if ismember(candid(j),model{1,I(i)})==1
                zz=find(model{1,I(i)}==candid(j));
                zz=zz(1);
                
                if zz>1
                    new_model{1,I(i)}=[model{1,I(i)}(1:zz-1),model{1,I(i)}(zz+1:end)];
                elseif length(model{1,I(i)})~=1
                    new_model{1,I(i)}=[model{1,I(i)}(2:end)];
                else
                    new_model{1,I(i)}=[];
                end
                
                %                 new_model{1,I(i)}=[model{1,I(i)}(1:zz-1),model{1,I(i)}(zz+1:end)];
                
                if length(new_model{1,I(i)})>=1
                    temp_var=ones(1,length(tf{1,1}));
                    
                    for jj=1:length(new_model{1,I(i)})
                        temp_var=temp_var.*tf{new_model{1,I(i)}(1,jj),1};
                    end
                else
                    temp_var=ones(1,length(tf{1,1}));
                end
                
                %
                %                 temp_var=ones(1,length(tf{1,1}));
                %                 alaki=0;
                %                 for jj=1:length(new_model{1,I(i)})
                %                     if new_model{1,I(i)}(1,jj)~=candid(j) | alaki==1
                %                         temp_var=temp_var.*tf{new_model{1,I(i)}(1,jj),1};
                %                     else
                %                         alaki=1;
                %                     end
                %                 end
                
                
                
                if temp_var~=1
                    new_variables{I(i),1}=temp_var;
                else
                    new_variables{I(i),1}=temp_var.*0;
                end
                [new_model,new_variables]=sing_check(new_model,new_variables);
%                 disp('model inside the model update22')
%                 
%                 
%                 for kkk=1:beta
%                     new_model{1,kkk}
%                 end
                [prid_error,coefs]=compute_error(new_variables,gene);
                error_energy(I(i),j+kk)=sum(sum(prid_error.^2));
            end
            
        end
    end
    
    % we have this error energy matrix to decide to jump
    
    error_energy;
    
    if min(min(error_energy))<c.*sum(sum(error{1,1}.^2))
        
        
        stop=0;
        [index1,index2]=find(error_energy==min(min(error_energy)));
        
        %         disp('before alaki')
        %
        %         error_energy
        %         min(min(error_energy))
        %         index1
        %         index2
        %
        %
        index1=index1(1);
        index2=index2(1);
        new_model=model;
        new_variables=variables;
        
        if index2<=kk
            new_model{1,index1}=[model{1,index1},candid(index2)];
            new_variables{index1,1}=new_variables{index1,1}.*tf{candid(index2),1};
        else
            %             disp('alaki')
            %             index1
            %             index2
            %             size(model)
            %             model{1,index1}
            %             candid(index2-kk)
            zz=find(model{1,index1}==candid(index2-kk));
            zz=zz(1);
            if zz>1
                new_model{1,index1}=[model{1,index1}(1:zz-1),model{1,index1}(zz+1:end)];
            elseif length(model{1,index1})~=1
                new_model{1,index1}=[model{1,index1}(2:end)];
            else
                new_model{1,index1}=[];
            end
            
            if length(new_model{1,index1})>=1
                temp_var=ones(1,length(tf{1,1}));
                
                for jj=1:length(new_model{1,index1})
                    temp_var=temp_var.*tf{new_model{1,index1}(1,jj),1};
                end
            else
                temp_var=ones(1,length(tf{1,1}));
            end
            
            if temp_var~=1
                new_variables{index1,1}=temp_var;
            else
                new_variables{index1,1}=temp_var.*0;
            end
        end
    else
        % the case where we don't have an improvement in error energy
        % we can try other jumps that we excluded in b for i=b+1:kk
        
        stop=1;
        new_model=model;
        new_variables=variables;
        
        
        
    end
    
end



