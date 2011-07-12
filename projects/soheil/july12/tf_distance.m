function d=tf_distance(tf1,tf2,ref)

len=length(ref);

index=find(ref==0);
loc1=index(1);
loc2=index(end);

[q1,a1,b1,m11,m21,m31]=seq_qual(tf1,ref);
[q2,a2,b2,m12,m22,m32]=seq_qual(tf2,ref);

% [loc1,loc2] determinis the middle interval

if loc1>1
%     disp('test')
%     sum(sum(tf1(1:loc1-1)~=tf2(1:loc1-1)))
%     if sum(sum(tf1(1:loc1-1)~=tf2(1:loc1-1)))~=0
%         [tf1(1:loc1-1)',tf2(1:loc1-1)']
%         temp=abs(corrcoef(tf1(1:loc1-1),tf2(1:loc1-1)))
%         d1=temp(1,2);
%     else
%         d1=1;
%     end
%     
    d1=sum(sum(tf1(1:loc1-1)==tf2(1:loc1-1)))/(loc1-1);
    
else
    disp('error in seq_qual')
end


% if sum(sum(tf1(loc1:loc2)~=tf2(loc1:loc2)))~=0
%     temp=abs(corrcoef(tf1(loc1:loc2),tf2(loc1:loc2)));
%     d2=temp(1,2);
% else
%     sum(sum(tf1(loc1:loc2)-tf2(loc1:loc2)))
%     d2=1;
% end

    d2=sum(sum(tf1(loc1:loc2)==tf2(loc1:loc2)))/(loc2-loc1+1);




if loc2<len
%     if sum(sum(tf1(loc2+1:end)~=tf2(loc2+1:end)))~=0
%         temp=abs(corrcoef(tf1(loc2+1:end),tf2(loc2+1:end)));
%         d3=temp(1,2);
%     else
%         d3=1;
%     end
    
        d3=sum(sum(tf1(loc2+1:end)==tf2(loc2+1:end)))/(len-loc2);
    
else
    disp('error in seq_qual')
end

if m11~=m12
   d1=0; 
end

if m21~=m22
   d2=0; 
end

if m31~=m32
   d3=0; 
end



d1
d2
d3
% d=(d1+d2+d3)/3

d=d1*d2*d3

