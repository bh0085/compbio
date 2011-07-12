function [q,aa,bb,m1,m2,m3,valid]=seq_qual(x,ref)

% both x and ref are trinary sequances of 0, 1 and -1
% ref is an increasing sequance (monotonic)

len=length(ref);

index=find(ref==0);
loc1=index(1);
loc2=index(end);
stop=0; % in a gene has a profile of two bits, need to be modified!
% [loc1,loc2] determinis the middle interval

if loc1>1
    [m1,f1]=mode(x(1:loc1-1));
else
    disp('error in seq_qual')
    loc1;
    stop=1;
end

[m2,f2]=mode(x(loc1:loc2));

if loc2<len
    [m3,f3]=mode(x(loc2+1:end));
else
    disp('error in seq_qual2')
    loc2
    stop=1;
end

% if seq is not monotonic, put quality=0;

if stop==0
    
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
    
    if q1<1/2
        q1=0;
    end
    
    if q2<1/2
        q2=0;
    end
    
    if q3<1/2
        q3=0;
    end
    
    
    %transition quality
%     if m1==m2
%         q12=0;
%     else
%         q12=1;
%     end
%     
%     if m2==m3
%         q23=0;
%     else
%         q23=1;
%     end
    
q12=1;
q23=1;
    
    aa=q1*q12*q2/2;
    bb=q2*q23*q3/2;
    
    q=(aa+bb)*check;
    valid=1;
else
    
    q=0;
    aa=0;
    bb=0;
    m1=0;
    m2=0;
    m3=0;
    valid=0;
    
end