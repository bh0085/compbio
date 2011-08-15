function [ind,pattern]=gene_pattern(ref)

% ref is an increasing sequance (monotonic)

len=length(ref);

index=find(ref==0);
loc1=index(1);
loc2=index(end);

if loc1<=floor(len/20) & loc2<len-floor(len/20)
    ind=1;
elseif loc2>=len-floor(len/20) & loc1>floor(len/20)
    ind=2;
elseif loc1<=floor(len/20) & loc2>=len-floor(len/20)
    ind=3;
else
    ind=0;
end

if ind==0
    pattern=[-1,0,1];
elseif ind==1
    pattern=[-1,0,-2];
elseif ind==2
    pattern=[0,1,-2];
else
%     disp('error in gene pattern')
    pattern=[-2,-2,-2];
end
