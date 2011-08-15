function [q,m1,m2,m3,ind]=seq_qual2(x,ref)

% both x and ref are trinary sequances of 0, 1 and -1
% ref is an increasing sequance (monotonic)

len=length(ref);

index=find(ref==0);
loc1=index(1);
loc2=index(end);

% if loc1<=floor(len/20) & loc2<len-floor(len/20)
%     ind=1;
% elseif loc2>=len-floor(len/20) & loc1>floor(len/20)
%     ind=2;
% elseif loc1<=floor(len/20) & loc2>=len-floor(len/20)
%     ind=3;
% else
%     ind=0;
% end
% 

[ind,pattern]=gene_pattern(ref);

if ind==0
    
    [m1,f1]=mode(x(1:loc1-1));
    
    [m2,f2]=mode(x(loc1:loc2));
    
    [m3,f3]=mode(x(loc2+1:end));
    
    % if the middle part is not 0, put quality=0; also exclude all 0 seq
    
    if m1==m3 & m1==0
        check=0;
    elseif m2~=0
        check=0;
    else
        check=1;
    end
    
    
    % state quality
    q1=f1/(loc1-1);
    q2=f2/(loc2-loc1+1);
    q3=f3/(len-loc2);
    
    % hard thresholding, at least majority should be more than half
    
    if q1<1/2
        q1=0;
    end
    
    if q2<1/2
        q2=0;
    end
    
    if q3<1/2
        q3=0;
    end
    
    
        q=(q1+q2+q3)*check/3;
    
%     q=(q1*q2+q2*q3)*check/2;
    
elseif ind==1
    
%     disp('length two pattern')
    
    [m1,f1]=mode(x(loc1:loc2));
    
    [m2,f2]=mode(x(loc2+1:end));
    
        m3=-2; % means it is not valid and the profile is just [m1,m2]

    
    % if the middle part is not 0, put quality=0; also exclude all 0 seq
    
    if m1==m2 & m1==0
        check=0;
    else
        check=1;
    end
    
    
    % state quality
    %     q1=f1/(loc1-1);
    q1=f1/(loc2-loc1+1);
    q2=f2/(len-loc2);
    
    % hard thresholding, at least majority should be more than half
    
    %     if q1<1/2
    %         q1=0;
    %     end
    
%     q1=-1; % not valid
    if q1<1/2
        q1=0;
    end
    
    if q2<1/2
        q2=0;
    end
    
    
        q=(q1+q2)*check/2;
    
        
%     q=(q2*q3)*check;
    
    
elseif ind==2
    
    
%     disp('length two pattern')
    [m1,f1]=mode(x(1:loc1-1));
    
    [m2,f2]=mode(x(loc1:loc2));
    
    m3=-2;
    
    % if the middle part is not 0, put quality=0; also exclude all 0 seq
    
    if m1==m2==0
        check=0;
    else
        check=1;
    end
    
    
    % state quality
    q1=f1/(loc1-1);
    q2=f2/(loc2-loc1+1);
    %     q3=f3/(len-loc2);
    
    % hard thresholding, at least majority should be more than half
    
    if q1<1/2
        q1=0;
    end
    
    if q2<1/2
        q2=0;
    end
    
    %     if q3<1/2
    %         q3=0;
    %     end
    %
    
    q3=-1;
        q=(q2+q1)*check/2;
    
%     q=(q2*q1)*check;
    
    
elseif ind==3
    
%     m1=-1;
    [m1,f1]=mode(x(loc1:loc2));
    m2=-2;
    m3=-2;
    
    % if the middle part is not 0, put quality=0; also exclude all 0 seq
    
    if m1==0
%         disp('constant profile')
    else
%         disp('error in gene quality profiling')
    end
    
    q1=-1;
    q2=-1;
    q3=-1;
    q=-1;
end


