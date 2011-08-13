function [var,m]=simplify_modules(var1,m1)

% this function simplifies the module codes
m1;
var=var1;
m=[];
index_or=find(m1==-1);

m_num=length(index_or)+1;

if length(m1)==0
var=var1;
m=m1;
elseif m_num==1 & length(m1)~=0
    
    %initialization
    m_temp=m1;
    
    len=length(m_temp);
    m_temp_new=[];
    if m_temp(1)~=-2
        m_temp_new(1)=m_temp(1);
        index=1;
    else
        m_temp_new(1:2)=m_temp(1:2);
        index=2;
    end
    not_point=0;
    
    %     disp('before')
    %     m_temp_new
    
    for j=index+1:len
        if m_temp(j)==-2
            not_point=1;
        else
            if not_point==0
                if ismember(m_temp(j),m_temp_new)==0
                    %                     disp('1')
                    m_temp_new=[m_temp_new,m_temp(j)];
                else
                    index_temp=find(m_temp_new==m_temp(j));
                    add=1;
                    
                    for ii=1:length(index_temp)
                        if index_temp(ii)==1
                            add=0;
                        elseif m_temp_new(index_temp(ii)-1)~=-2
                            add=0;
                        end
                    end
                    
                    if add==1
                        m_temp_new=[m_temp_new,m_temp(j)];
                    end
                    %                         index_temp=index_temp(1);
                    %                         if index_temp>=2 & m_temp_new(index_temp-1)==-2
                    %                             disp('2')
                    %                             m_temp_new=[m_temp_new,m_temp(j)]
                    %                         end
                    
                    
                end
            else
                if ismember(m_temp(j),m_temp_new)==0
                    %                     disp('3')
                    m_temp_new=[m_temp_new,-2,m_temp(j)];
                else
                    index_temp=find(m_temp_new==m_temp(j));
                    
                    add=1;
                    
                    for ii=1:length(index_temp)
                        
                        if index_temp(ii)>=2 & m_temp_new(index_temp(ii)-1)==-2
                            add=0;
                        end
                    end
                    
                    if add==1
                        m_temp_new=[m_temp_new,-2,m_temp(j)];
                    end
                    
                    %                     index_temp=index_temp(1);
                    %                     if index_temp>=2 & m_temp_new(index_temp-1)==-2
                    %                         % dont add
                    %                     else
                    %                         disp('4')
                    %                         m_temp_new=[m_temp_new,-2,m_temp(j)]
                    %                     end
                end
            end
            not_point=0;
            
        end
        
    end
    
%     disp('here')
    m=m_temp_new;
    
else
    
    for i=1:m_num
        if i==1
            m_temp=m1(1:index_or(1)-1);
        elseif i==m_num
            m_temp=m1(index_or(end)+1:end);
        else
            m_temp=m1(index_or(i-1)+1:index_or(i)-1);
        end
        
        m_temp;
        
        len=length(m_temp);
        if len~=0
        m_temp_new=[];
        if m_temp(1)~=-2
            m_temp_new(1)=m_temp(1);
            index=1;
        else
            m_temp_new(1:2)=m_temp(1:2);
            index=2;
        end
        not_point=0;
        
        for j=index+1:len
            if m_temp(j)==-2
                not_point=1;
            else
                if not_point==0
                    if ismember(m_temp(j),m_temp_new)==0
                        m_temp_new=[m_temp_new,m_temp(j)];
                    else
                        
                        index_temp=find(m_temp_new==m_temp(j));
                        %                         index_temp=index_temp(1);
                        
                        add=1;
                        
                        for ii=1:length(index_temp)
                            if index_temp(ii)==1
                                add=0;
                            elseif m_temp_new(index_temp(ii)-1)~=-2
                                add=0;
                            end
                        end
                        
                        if add==1
                            m_temp_new=[m_temp_new,m_temp(j)];
                        end
                        
                        %                         if index_temp>=2 & m_temp_new(index_temp-1)==-2
                        %
                        %                             m_temp_new=[m_temp_new,m_temp(j)];
                        %                         end
                    end
                else
                    if ismember(m_temp(j),m_temp_new)==0
                        m_temp_new=[m_temp_new,-2,m_temp(j)];
                    else
                        index_temp=find(m_temp_new==m_temp(j));
                        
                        add=1;
                        
                        for ii=1:length(index_temp)
                            
                            if index_temp(ii)>=2 & m_temp_new(index_temp(ii)-1)==-2
                                add=0;
                            end
                        end
                        
                        if add==1
                            m_temp_new=[m_temp_new,-2,m_temp(j)];
                        end
                        %                         index_temp=index_temp(1);
                        %                         if index_temp>=2 & m_temp_new(index_temp-1)==-2
                        %                             % dont add
                        %                         else
                        %                             m_temp_new=[m_temp_new,-2,m_temp(j)];
                        %                         end
                    end
                end
                not_point=0;
                
            end
            
        end
        m=[m,m_temp_new,-1];
    end
    end
%     disp('2')
    m=m(1:end-1);
end

% elimination same or modules
index_or=find(m==-1);
num_m=length(find(m==-1))+1;


if num_m>=2
    mm=cell(1,1);
    mm{1,1}=[m(1:index_or(1)-1)];
    for i=2:num_m
        num_temp=size(mm,1);
        if i~=num_m
            cand_mod=m(index_or(i-1)+1:index_or(i)-1);
        else
            cand_mod=m(index_or(i-1)+1:end);
        end
        
        add=1;
        for j=1:num_temp
            if length(mm{j,1})==length(cand_mod)
                if length(intersect(mm{j,1},cand_mod))==length(cand_mod)
                    add=0;
                end
            end
        end
        
        if add==1
            mm{num_temp+1,1}=cand_mod;
        end
    end
    m=[];
    for i=1:size(mm,1)
        m=[m,mm{i,1},-1];
    end
    m=m(1:end-1);
end













