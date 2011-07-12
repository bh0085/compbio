function [var,m]=my_and_modules(var1,m1,var2,m2)

var=union(var1,var2);

index1=find(m1==-1);
index2=find(m2==-1);

m1_num=length(index1)+1;
m2_num=length(index2)+1;

m_num=m1_num*m2_num;
m_cell=cell(1,m_num); % each cell has between two -1 (like monomials)


for j1=1:m1_num
    for j2=1:m2_num
        if m_num==1
            m_cell{1,(j1-1)*m2_num+j2}= [m1,m2];
            
        elseif m1_num==1
            
            if j2==1
                m2_temp=m2(1:index2(1)-1);
            elseif j2==m2_num
                m2_temp=m2(index2(end)+1:end);
            else
                m2_temp=m2(index2(j2-1)+1:index2(j2)-1);
            end
            
            m_cell{1,(j1-1)*m2_num+j2}=[m1,m2_temp];
            
        elseif m2_num==1
            
            
             if j1==1
                m1_temp=m1(1:index1(1)-1);
            elseif j1==m1_num
                m1_temp=m1(index1(end)+1:end);
            else
                m1_temp=m1(index1(j1-1)+1:index1(j1)-1);
            end
            
            m_cell{1,(j1-1)*m2_num+j2}=[m1_temp,m2];
            
        else
            
             if j1==1
                m1_temp=m1(1:index1(1)-1);
            elseif j1==m1_num
                m1_temp=m1(index1(end)+1:end);
            else
                m1_temp=m1(index1(j1-1)+1:index1(j1)-1);
            end
            
            
            if j2==1
                m2_temp=m2(1:index2(1)-1);
            elseif j2==m2_num
                m2_temp=m2(index2(end)+1:end);
            else
                m2_temp=m2(index2(j2-1)+1:index2(j2)-1);
            end
                        
           m_cell{1,(j1-1)*m2_num+j2}=[m1_temp,m2_temp];

        end
        
        
    end
end


m=[];
for i=1:m_num
    
   m=[m,m_cell{1,i},-1]; 
end

m=m(1,1:end-1);

