function [index1,index2]=rand_walk(temp)

% disp('check validity of the stoch matrix')
% sum(sum(temp))


r=rand(1,1);
[k,num_tf]=size(temp);
s1=0;
s2=0;
for i=1:k
    for j=1:num_tf
        
        s1=s2;
        s2=s1+temp(i,j);
%         
%         i
%         j
%         s1
%         s2
%         temp(i,j)
        
        if r>s1 & r<s2
%             disp('inside rand walk')
%             i
%             j
%             s1
%             s2

            index1=i;
            index2=j;
        end
    end
end