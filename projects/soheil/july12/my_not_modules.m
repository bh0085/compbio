function [var,m_new]=my_not_modules(var,m)

if length(m)~=0

index=find(m==-1);

m_num=length(index)+1;

m_cell=cell(1,m_num);

for j=1:m_num
    
    if m_num==1
        m_cell{1,1}=[];
        not_flag=0;
        for i=1:length(m)
            if m(i)>0 & not_flag==0
                m_cell{1,j}=[m_cell{1,j},-2,m(i),-1];
                not_flag=0;
            elseif m(i)>0 & not_flag==1
                m_cell{1,j}=[m_cell{1,j},m(i),-1];
                not_flag=0;
            else
                not_flag=1;
            end
%             disp('inside not')
%             not_flag
%             m_cell{1,j}
        end
                m_cell{1,j}=m_cell{1,j}(1,1:end-1);

        
    elseif j==1
        m_cell{1,1}=[];
        not_flag=0;
        start_ind=1;
        end_ind=index(1)-1;
        
        for i=start_ind:end_ind
            if m(i)>0 & not_flag==0
                m_cell{1,j}=[m_cell{1,j},-2,m(i),-1];
                not_flag=0;
            elseif m(i)>0 & not_flag==1
                m_cell{1,j}=[m_cell{1,j},m(i),-1];
                not_flag=0;
            else
                not_flag=1;
            end
        end
        
        m_cell{1,j}=m_cell{1,j}(1,1:end-1);
        
    elseif j==m_num
        
        m_cell{1,j}=[];
        not_flag=0;
        start_ind=index(end)+1;
        end_ind=length(m);
        
        for i=start_ind:end_ind
            if m(i)>0 & not_flag==0
                m_cell{1,j}=[m_cell{1,j},-2,m(i),-1];
                not_flag=0;
            elseif m(i)>0 & not_flag==1
                m_cell{1,j}=[m_cell{1,j},m(i),-1];
                not_flag=0;
            else
                not_flag=1;
            end
        end
        
        m_cell{1,j}=m_cell{1,j}(1,1:end-1);

        
    else
        
        m_cell{1,j}=[];
        not_flag=0;
        start_ind=index(j-1)+1;
        end_ind=index(j)-1;
        
        for i=start_ind:end_ind
            if m(i)>0 & not_flag==0
                m_cell{1,j}=[m_cell{1,j},-2,m(i),-1];
                not_flag=0;
            elseif m(i)>0 & not_flag==1
                m_cell{1,j}=[m_cell{1,j},m(i),-1];
                not_flag=0;
            else
                not_flag=1;
            end
        end
        
              m_cell{1,j}=m_cell{1,j}(1,1:end-1);
  
    end
    
    
end

m_new=[];
if m_num==1
    m_new=m_cell{1,1};
else
    for j=1:m_num
        [temp,m_new]=my_and_modules(var,m_new,var,m_cell{1,j});
    end
end

else
    m_new=m;
end

