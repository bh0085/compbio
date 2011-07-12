function   find_logic_model(gene,tf_sq,tf_list,tf_and_list,tf_or_list)

th=0.6;
pf_d01=[];
pf_001=[];
pf_10d=[];
pf_100=[];
pf_00d=[];
pf_d00=[];

for i=1:size(tf_list,2)
    
    if tf_list(5:7,i)==[-1 0 1]'
        pf_d01=[pf_d01,[tf_list(1,i),0,0,tf_list(2:4,i)']'];
        % first two rows are TF names
        % third row is a control bit: 0: ind TF/ 2: or/1: and, etc
        % the next three rows are qualities
        
    elseif tf_list(5:7,i)==[0 0 1]'
        pf_001=[pf_001,[tf_list(1,i),0,0,tf_list(2:4,i)']'];
    elseif tf_list(5:7,i)==[1 0 -1]'
        pf_10d=[pf_10d,[tf_list(1,i),0,0,tf_list(2:4,i)']'];
    elseif tf_list(5:7,i)==[1 0 0]'
        pf_100=[pf_100,[tf_list(1,i),0,0,tf_list(2:4,i)']'];
    elseif tf_list(5:7,i)==[0 0 -1]'
        pf_00d=[pf_00d,[tf_list(1,i),0,0,tf_list(2:4,i)']'];
    elseif tf_list(5:7,i)==[-1 0 0]'
        pf_d00=[pf_d00,[tf_list(1,i),0,0,tf_list(2:4,i)']'];
    else
        disp('unknown profile')
        tf_list(5:7,i)
    end
    
end


% updating profiles by and modules

for i=1:size(tf_and_list,2)
    
    if tf_and_list(6:8,i)==[-1 0 1]'
        pf_d01=[pf_d01,[tf_and_list(1,i),tf_and_list(2,i),1,tf_and_list(3:5,i)']'];
        
    elseif tf_and_list(6:8,i)==[0 0 1]'
        pf_001=[pf_001,[tf_and_list(1,i),tf_and_list(2,i),1,tf_and_list(3:5,i)']'];
    elseif tf_and_list(6:8,i)==[1 0 -1]'
        pf_10d=[pf_10d,[tf_and_list(1,i),tf_and_list(2,i),1,tf_and_list(3:5,i)']'];
    elseif tf_and_list(6:8,i)==[1 0 0]'
        pf_100=[pf_100,[tf_and_list(1,i),tf_and_list(2,i),1,tf_and_list(3:5,i)']'];
    elseif tf_and_list(6:8,i)==[0 0 -1]'
        pf_00d=[pf_00d,[tf_and_list(1,i),tf_and_list(2,i),1,tf_and_list(3:5,i)']'];
    elseif tf_and_list(6:8,i)==[-1 0 0]'
        pf_d00=[pf_d00,[tf_and_list(1,i),tf_and_list(2,i),1,tf_and_list(3:5,i)']'];
    else
        disp('unknown profile')
        tf_and_list(6:8,i)
    end
    
end


for i=1:size(tf_or_list,2)
    
    if tf_or_list(6:8,i)==[-1 0 1]'
        pf_d01=[pf_d01,[tf_or_list(1,i),tf_or_list(2,i),1,tf_or_list(3:5,i)']'];
        
    elseif tf_or_list(6:8,i)==[0 0 1]'
        pf_001=[pf_001,[tf_or_list(1,i),tf_or_list(2,i),2,tf_or_list(3:5,i)']'];
    elseif tf_or_list(6:8,i)==[1 0 -1]'
        pf_10d=[pf_10d,[tf_or_list(1,i),tf_or_list(2,i),2,tf_or_list(3:5,i)']'];
    elseif tf_or_list(6:8,i)==[1 0 0]'
        pf_100=[pf_100,[tf_or_list(1,i),tf_or_list(2,i),2,tf_or_list(3:5,i)']'];
    elseif tf_or_list(6:8,i)==[0 0 -1]'
        pf_00d=[pf_00d,[tf_or_list(1,i),tf_or_list(2,i),2,tf_or_list(3:5,i)']'];
    elseif tf_or_list(6:8,i)==[-1 0 0]'
        pf_d00=[pf_d00,[tf_or_list(1,i),tf_or_list(2,i),2,tf_or_list(3:5,i)']'];
    else
        disp('unknown profile')
        tf_or_list(6:8,i)
    end
    
end

if size(pf_10d,2)~=0
    pf_10d(3,:)=pf_10d(3,:)+3;
end
% 3:not ind
% 4: not and
% 5: not or

X11=[pf_d01,pf_10d];
X1=X11;
if size(X1,2)~=0
    [temp,I]=sort(X11(4,:),'descend');
    for j=1:size(X1,2)
        X1(:,j)=X11(:,I(j));
    end
    
    index=find(X1(4,:)>th);
    if length(index)~=0
        disp('X1')
        ind_X1=index(end);
        
        X1(:,1:ind_X1);
    end
    
end




if size(pf_00d,2)~=0
    pf_00d(3,:)=pf_00d(3,:)+3;
end

X22=[pf_001,pf_00d];



X2=X22;
if size(X2,2)~=0
    [temp,I]=sort(X22(4,:),'descend');
    for j=1:size(X2,2)
        X2(:,j)=X22(:,I(j));
    end
end



if size(pf_100,2)~=0
    pf_100(3,:)=pf_100(3,:)+3;
end
X33=[pf_d00,pf_100];

X3=X33;
if size(X3,2)~=0
    [temp,I]=sort(X33(4,:),'descend');
    for j=1:size(X3,2)
        X3(:,j)=X33(:,I(j));
    end
end


% X2
% X3

% X12:
% X12: X1 and X2
% Y12: X1 or X2
% first 6 rows are copies from the first three rows of X1 and X2
%the next three rows are qualities
len12=size(X1,2)*size(X2,2);
if len12>0
    X12=zeros(9,len12);
    Y12=zeros(9,len12);
    
    
    for i=1:size(X1,2)
        for j=1:size(X2,2)
            X12(1:3,(i-1)*size(X1,2)+j)=X1(1:3,i);
            X12(4:6,(i-1)*size(X1,2)+j)=X2(1:3,j);
            
            Y12(1:3,(i-1)*size(X1,2)+j)=X1(1:3,i);
            Y12(4:6,(i-1)*size(X1,2)+j)=X2(1:3,j);
            
            if X1(3,i)==0 & X2(3,j)==0
                [q,q1,q2,p1,p2,p3,valid]=seq_qual(my_and(tf_sq(X1(1,i),:),tf_sq(X2(1,j),:)),gene);
                X12(7:9,(i-1)*size(X1,2)+j)=[q,q1,q2]';
                
                [q,q1,q2,p1,p2,p3,valid]=seq_qual(my_or(tf_sq(X1(1,i),:),tf_sq(X2(1,j),:)),gene);
                Y12(7:9,(i-1)*size(X1,2)+j)=[q,q1,q2]';
                
                
            elseif X1(3,i)==0 & X2(3,j)==3
                [q,q1,q2,p1,p2,p3,valid]=seq_qual(my_and(tf_sq(X1(1,i),:),my_not(tf_sq(X2(1,j),:))),gene);
                X12(7:9,(i-1)*size(X1,2)+j)=[q,q1,q2]';
                
                [q,q1,q2,p1,p2,p3,valid]=seq_qual(my_or(tf_sq(X1(1,i),:),my_not(tf_sq(X2(1,j),:))),gene);
                Y12(7:9,(i-1)*size(X1,2)+j)=[q,q1,q2]';
                
            elseif X1(3,i)==3 & X2(3,j)==0
                [q,q1,q2,p1,p2,p3,valid]=seq_qual(my_and(my_not(tf_sq(X1(1,i),:)),tf_sq(X2(1,j),:)),gene);
                X12(7:9,(i-1)*size(X1,2)+j)=[q,q1,q2]';
                
                [q,q1,q2,p1,p2,p3,valid]=seq_qual(my_or(my_not(tf_sq(X1(1,i),:)),tf_sq(X2(1,j),:)),gene);
                Y12(7:9,(i-1)*size(X1,2)+j)=[q,q1,q2]';
                
            elseif X1(3,i)==3 & X2(3,j)==3
                [q,q1,q2,p1,p2,p3,valid]=seq_qual(my_and(my_not(tf_sq(X1(1,i),:)),my_not(tf_sq(X2(1,j),:))),gene);
                X12(7:9,(i-1)*size(X1,2)+j)=[q,q1,q2]';
                
                [q,q1,q2,p1,p2,p3,valid]=seq_qual(my_or(my_not(tf_sq(X1(1,i),:)),my_not(tf_sq(X2(1,j),:))),gene);
                Y12(7:9,(i-1)*size(X1,2)+j)=[q,q1,q2]';
                
                %             elseif X1(3,i)==1 & X2(3,j)==1
                %                 [q,q1,q2,p1,p2,p3,valid]=seq_qual(my_and(tf_sq(X1(1,i),:),tf_sq(X2(1,j),:)),gene);
                %                 X12(7:9,(i-1)*size(X1,2)+j)=[q,q1,q2]';
                %
                %             elseif X1(3,i)==1 & X2(3,j)==2
                %                 [q,q1,q2,p1,p2,p3,valid]=seq_qual(my_and(tf_sq(X1(1,i),:),tf_sq(X2(1,j),:)),gene);
                %                 X12(7:9,(i-1)*size(X1,2)+j)=[q,q1,q2]';
                %
                %             elseif X1(3,i)==2 & X2(3,j)==0
                %                 [q,q1,q2,p1,p2,p3,valid]=seq_qual(my_and(tf_sq(X1(1,i),:),tf_sq(X2(1,j),:)),gene);
                %                 X12(7:9,(i-1)*size(X1,2)+j)=[q,q1,q2]';
                %
                %             elseif X1(3,i)==2 & X2(3,j)==1
                %                 [q,q1,q2,p1,p2,p3,valid]=seq_qual(my_and(tf_sq(X1(1,i),:),tf_sq(X2(1,j),:)),gene);
                %                 X12(7:9,(i-1)*size(X1,2)+j)=[q,q1,q2]';
                %
                %             elseif X1(3,i)==2 & X2(3,j)==2
                %                 [q,q1,q2,p1,p2,p3,valid]=seq_qual(my_and(tf_sq(X1(1,i),:),tf_sq(X2(1,j),:)),gene);
                %                 X12(7:9,(i-1)*size(X1,2)+j)=[q,q1,q2]';
                
                
            end
        end
        
    end
    
    
    X12_s=X12;
    [temp,I]=sort(X12(7,:),'descend');
    for j=1:size(X12,2)
        X12_s(:,j)=X12(:,I(j));
    end
    
    X12=X12_s;
    
    
    Y12_s=Y12;
    [temp,I]=sort(Y12(7,:),'descend');
    for j=1:size(Y12,2)
        Y12_s(:,j)=Y12(:,I(j));
    end
    
    Y12=Y12_s;
    
    
    index=find(X12(7,:)>th);
    if length(index)~=0
        disp('X12')
        ind_X12=index(end);
        X12(:,1:ind_X12)
    end
    
    index=find(Y12(7,:)>th);
    
    if length(index)~=0
        disp('Y12')
        ind_Y12=index(end);
        
        Y12(:,1:ind_Y12)
    end
end



% X13: X1 and X3
% Y13: Y1 or Y3

len13=size(X1,2)*size(X3,2);
if len13>0
    X13=zeros(9,len13);
    Y13=zeros(9,len13);
    
    
    for i=1:size(X1,2)
        for j=1:size(X3,2)
            X13(1:3,(i-1)*size(X1,2)+j)=X1(1:3,i);
            X13(4:6,(i-1)*size(X1,2)+j)=X3(1:3,j);
            
            Y13(1:3,(i-1)*size(X1,2)+j)=X1(1:3,i);
            Y13(4:6,(i-1)*size(X1,2)+j)=X3(1:3,j);
            
            if X1(3,i)==0 & X3(3,j)==0
                [q,q1,q2,p1,p2,p3,valid]=seq_qual(my_and(tf_sq(X1(1,i),:),tf_sq(X3(1,j),:)),gene);
                X13(7:9,(i-1)*size(X1,2)+j)=[q,q1,q2]';
                
                [q,q1,q2,p1,p2,p3,valid]=seq_qual(my_or(tf_sq(X1(1,i),:),tf_sq(X3(1,j),:)),gene);
                Y13(7:9,(i-1)*size(X1,2)+j)=[q,q1,q2]';
                
            elseif X1(3,i)==0 & X3(3,j)==3
                [q,q1,q2,p1,p2,p3,valid]=seq_qual(my_and(tf_sq(X1(1,i),:),my_not(tf_sq(X3(1,j),:))),gene);
                X13(7:9,(i-1)*size(X1,2)+j)=[q,q1,q2]';
                
                [q,q1,q2,p1,p2,p3,valid]=seq_qual(my_or(tf_sq(X1(1,i),:),my_not(tf_sq(X3(1,j),:))),gene);
                Y13(7:9,(i-1)*size(X1,2)+j)=[q,q1,q2]';
                
            elseif X1(3,i)==3 & X3(3,j)==0
                [q,q1,q2,p1,p2,p3,valid]=seq_qual(my_and(my_not(tf_sq(X1(1,i),:)),tf_sq(X3(1,j),:)),gene);
                X13(7:9,(i-1)*size(X1,2)+j)=[q,q1,q2]';
                
                [q,q1,q2,p1,p2,p3,valid]=seq_qual(my_or(my_not(tf_sq(X1(1,i),:)),tf_sq(X3(1,j),:)),gene);
                Y13(7:9,(i-1)*size(X1,2)+j)=[q,q1,q2]';
                
            elseif X1(3,i)==3 & X3(3,j)==3
                [q,q1,q2,p1,p2,p3,valid]=seq_qual(my_and(my_not(tf_sq(X1(1,i),:)),my_not(tf_sq(X3(1,j),:))),gene);
                X13(7:9,(i-1)*size(X1,2)+j)=[q,q1,q2]';
                
                [q,q1,q2,p1,p2,p3,valid]=seq_qual(my_or(my_not(tf_sq(X1(1,i),:)),my_not(tf_sq(X3(1,j),:))),gene);
                Y13(7:9,(i-1)*size(X1,2)+j)=[q,q1,q2]';
                
            end
        end
        
    end
    
    
    X13_s=X13;
    [temp,I]=sort(X13(7,:),'descend');
    for j=1:size(X13,2)
        X13_s(:,j)=X13(:,I(j));
    end
    
    X13=X13_s;
    
    Y13_s=Y13;
    [temp,I]=sort(Y13(7,:),'descend');
    for j=1:size(Y13,2)
        Y13_s(:,j)=Y13(:,I(j));
    end
    
    
    Y13=Y13_s;
    
    
    
    index=find(X13(7,:)>th);
    if length(index)~=0
        disp('X13')
        ind_X13=index(end);
        X13(:,1:ind_X13)
    end
    
    index=find(Y13(7,:)>th);
    if length(index)~=0
        disp('Y13')
        ind_Y13=index(end);
        
        Y13(:,1:ind_Y13)
    end
end



% Y23=Y2+Y3

len23=size(X2,2)*size(X3,2);
if len23>0
    Y23=zeros(9,len23);
    
    
    for i=1:size(X2,2)
        for j=1:size(X3,2)
            Y23(1:3,(i-1)*size(X2,2)+j)=X2(1:3,i);
            Y23(4:6,(i-1)*size(X2,2)+j)=X3(1:3,j);
            
            if X2(3,i)==0 & X3(3,j)==0
                [q,q1,q2,p1,p2,p3,valid]=seq_qual(my_or(tf_sq(X2(1,i),:),tf_sq(X3(1,j),:)),gene);
                Y23(7:9,(i-1)*size(X2,2)+j)=[q,q1,q2]';
                
            elseif X2(3,i)==0 & X3(3,j)==3
                [q,q1,q2,p1,p2,p3,valid]=seq_qual(my_or(tf_sq(X2(1,i),:),my_not(tf_sq(X3(1,j),:))),gene);
                Y23(7:9,(i-1)*size(X2,2)+j)=[q,q1,q2]';
                
            elseif X2(3,i)==3 & X3(3,j)==0
                [q,q1,q2,p1,p2,p3,valid]=seq_qual(my_or(my_not(tf_sq(X2(1,i),:)),tf_sq(X3(1,j),:)),gene);
                Y23(7:9,(i-1)*size(X2,2)+j)=[q,q1,q2]';
                
            elseif X2(3,i)==3 & X3(3,j)==3
                [q,q1,q2,p1,p2,p3,valid]=seq_qual(my_or(my_not(tf_sq(X2(1,i),:)),my_not(tf_sq(X3(1,j),:))),gene);
                Y23(7:9,(i-1)*size(X2,2)+j)=[q,q1,q2]';
                
            end
        end
        
    end
    
    
    Y23_s=Y23;
    [temp,I]=sort(Y23(7,:),'descend');
    for j=1:size(Y23,2)
        Y23_s(:,j)=Y23(:,I(j));
    end
    
    Y23=Y23_s;
    
    index=find(Y23(7,:)>th);
    
    if length(index)~=0
        
        ind_Y23=index(end);
        disp('Y23')
        Y23(:,1:ind_Y23)
    end
end





