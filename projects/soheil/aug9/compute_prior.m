function p=compute_prior(prior,module)

% this function computes complexity penalty of a circuit
% prior is a scalar (say 0.75) penalizing each gate
% module is the code for that circuit

index=find(module<0);
gate_num=length(index);

index_or=find(module==-1);
num_module=length(index_or)+1;

if num_module==1
    m_temp=module;
    
    % compute the prior
    num_not=length(find(m_temp==-2));
    gate_num=gate_num+length(m_temp)-num_not-1;
    
else
    for i=1:num_module
        if i==1
            m_temp=module(1:index_or(1)-1);
        elseif i==num_module
            m_temp=module(index_or(end)+1:end);
        else
            m_temp=module(index_or(i-1)+1:index_or(i)-1);
        end
        
        % compute the prior
        num_not=length(find(m_temp==-2));
        gate_num=gate_num+length(m_temp)-num_not-1;
        
    end
    
end

gate_num;

p=prior^gate_num;
