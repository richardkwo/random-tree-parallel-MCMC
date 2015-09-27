function averaged_chain = aggregate_average(sub_chain)
fprintf('Combining chains with simple averaging...\n');
% aggregation via simple averaging
k = length(sub_chain);
[~, p] = size(sub_chain{1});
chain_lengths = zeros(1,k);
for c=1:k
    chain_lengths(c) = size(sub_chain{c},1);
end
n = min(chain_lengths);
averaged_chain = zeros(n, p);
for c=1:k
    averaged_chain = averaged_chain + sub_chain{c}(1:n, :) * 1/k;
end

end