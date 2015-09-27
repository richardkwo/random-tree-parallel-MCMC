function [index, value, goodness_of_cuts] = consensusMLCut(x, l, r, M, subset_index)
% [index, value] = consensusMLCut(x, l, r, M, subset_index)
% get the cut that maximizes the summed empirical likelihood from multiple subsets via a line search
% 
% ARGUMENT
% x: 1-D array
% l: scalar, left boundary
% r: scalar, right boundary
% M: scalar, # of subsets
% subset_index: (i.e. mark(:,2)) the indicator for which subset the sample comes from, of the same 
% dimension as x
% 
% OUTPUT
% index: the index of x that is the rightmost point lying left of cutting point
% value: the value of the cutting point
%
x = x(:)';
N = length(x);
[y, I] = sort(x); % increasing order
assert(length(x) == length(subset_index), 'length of x must match length of subset_index');
assert(all(subset_index>=1) && all(subset_index<=M), 'subset index must take 1...M');
assert(y(1)>=l & y(end)<=r, 'x not in range of [l,r]');
% index_vec is a subset of I and is considered for cutting
index_vec = 1:(N-1); % corresponding to x(index_vec(1:end-1))
cuts = y(index_vec);
num_cuts = length(index_vec);
subset_sizes = histc(subset_index, 1:M);
% profile a (M x num_cuts) cum-sum matrix: cum_sum_matrix(m, j) = # of samples in subset m that falls left of cut j
sorted_subset_indices = subset_index(I(index_vec));
cum_sum_matrix = zeros(M, num_cuts);
complement_cum_sum_matrix = zeros(M, num_cuts);
for m=1:M
    cum_sum_matrix(m,:) = cumsum(sorted_subset_indices==m); % deal with zeros 
    complement_cum_sum_matrix(m,:) = subset_sizes(m) - cum_sum_matrix(m,:);
end
% compute goodness of cut from matrix operations
L = cuts - l;
R = r - cuts;
L(L==0) = 1e-4;
R(R==0) = 1e-4;
A1 = repmat(L, M, 1);
A2 = repmat(R, M, 1);
tmp_mat_1 = cum_sum_matrix .* (log(cum_sum_matrix) - log(A1));
tmp_mat_1(cum_sum_matrix==0) = 0;
tmp_mat_2 = complement_cum_sum_matrix .* (log(complement_cum_sum_matrix) - log(A2));
tmp_mat_2(complement_cum_sum_matrix==0) = 0;
goodness_of_cuts = sum(tmp_mat_1 + tmp_mat_2, 1);
% pick the ML
[~, u] = max(goodness_of_cuts);
index = I(u);
value = x(index);
% fprintf('%d: goc = %f\n', u, goodness_of_cuts(u));
end