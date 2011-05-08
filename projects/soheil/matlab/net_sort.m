function [net_sorted,net_prob_sorted]=net_sort(net_prob)

[g_num,tf_num]=size(net_prob);
net_sorted=net_prob.*0;
net_prob_sorted=net_prob.*0;

for i=1:g_num
    [sorted_row,I] = sort(net_prob(i,:),'descend');
    net_sorted(i,:)=I;
    net_prob_sorted(i,:)=sorted_row;
end
