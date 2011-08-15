function y=read_module(tf_id,code_id,tf_sq)

% tf_if is a vector that says the id of TF's e.g.: [2 5 41]
% code id: uses TF names
% for not of a variable: uses -3 before it, not only applies on single variables
% or is -1 , NO AND SIGN, between two -1 everything is anded
%example: TF4+TF4*not(TF2)+TF2*TF3 =>
% [TF4 TF2 TF3], [4 -1 4 -3 2 -1 2 3]
% tf_sq is a matrix of size tf_num* len
% y is a vector of size 1* len


[tf_num,len]=size(tf_sq);
y=zeros(1,len);

c_len=length(code_id);
var_num=length(tf_id);

index=find(code_id==-1);
y_modules=zeros(length(index)+1,len);

if length(index)==0
    
    len_temp=c_len;
    not_sign=0;
    for i=1:len_temp
        if code_id(i)==-2
            not_sign=1;
        elseif not_sign==0

            not_sign=0;
            y_modules(1,:)=my_and(y_modules(1,:),tf_sq(code_id(i),:));
            
        else
            not_sign=0;
            y_modules(1,:)=my_and(y_modules(1,:),my_not(tf_sq(code_id(i),:)));
        end
        
    end
    
else
    
    for j=1:length(index)+1
        
        if j>=2 & j~=length(index)+1
            len_temp=index(j)-index(j-1)-1;
            start_index=index(j-1)+1;
        elseif j>=2
            len_temp=c_len-index(j-1);
            start_index=index(j-1)+1;
            
        else
            len_temp=index(j)-1;
            start_index=1;
        end
        
        not_sign=0;
        for i=1:len_temp
            if code_id(start_index+i-1)==-2
                not_sign=1;
            elseif not_sign==0
                not_sign=0;
                y_modules(j,:)=my_and(y_modules(j,:),tf_sq(code_id(start_index+i-1),:));
            else
                not_sign=0;
                y_modules(j,:)=my_and(y_modules(j,:),my_not(tf_sq(code_id(start_index+i-1),:)));
            end
            
        end
        
        
    end
    
    
    
end



y=zeros(1,len);
for j=1:length(index)+1
    y=my_or(y,y_modules(j,:));
end








