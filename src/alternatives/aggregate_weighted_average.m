function averaged_chain = aggregate_weighted_average(sub_chain)
fprintf('Aggregating by weighted averaging...\n');
k = length(sub_chain);
p = size(sub_chain{1},2);
chain_length = zeros(1,k);
for c=1:k
    chain_length(c) = size(sub_chain{c},1);
end
n = min(chain_length);
averaged_chain = zeros(n, p);

w = cell(1,k);
w_inv = cell(1,k);
sum_inv = zeros(p);
for i = 1:k
    w{i} = cov(sub_chain{i});
    w_inv{i} = inv(w{i});
    sum_inv = sum_inv + w_inv{i};
end

for c=1:k
    averaged_chain = averaged_chain + sub_chain{c}(1:n,:) * w_inv{i} / sum_inv;
end


end