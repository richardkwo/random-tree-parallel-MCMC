function aggregated_samples = aggregate_uai_nonparametric(sub_chain, varargin)
% drawing N posterior samples from multiplicated KD estimated subchain posteriors
% "nonparametric" method from the uai paper. See its Algorithm 1. 

if ~isempty(varargin)
    N = varargin{1};
else
    N = size(sub_chain{1},1);
end

M = length(sub_chain);
p = size(sub_chain{1},2);
subchain_sizes = zeros(1,M);
for c=1:M
    subchain_sizes(c) = size(sub_chain{c},1);
end

aggregated_samples = zeros(N, p);
index_t = zeros(1, M);
% uniformly init the index 
for c=1:M
    index_t(c) = randsample(subchain_sizes(c), 1);
end
x_t = zeros(M, p);
for c=1:M
    x_t(c,:) = sub_chain{c}(index_t(c), :);
end
x_t_mean = mean(x_t, 1);
ss_t = sum(sum((x_t - repmat(x_t_mean, M, 1)).^2));

fprintf('Getting %d samples with nonparametric method...\n', N);
tic;
count_rej = 0;
for i=1:N
    h = i^(-1/(4+p));
    for c=1:M
        index_new = index_t;
        % move the index for this subchain
        index_new(c) = randsample(subchain_sizes(c), 1); 
        x_t_new = x_t;
        x_t_new(c,:) = sub_chain{c}(index_new(c), :);
        x_t_new_mean = mean(x_t_new, 1);
        ss_new = sum(sum((x_t_new - repmat(x_t_new_mean, M, 1)).^2));
        u = rand;
        if u < exp(-(ss_new - ss_t)/(2*h^2))
            % accept the move
            index_t = index_new;
            x_t = x_t_new;
            x_t_mean = x_t_new_mean;
            ss_t = ss_new;
        else
            count_rej = count_rej + 1;
        end
    end
    % sample with the index
    aggregated_samples(i, :) = mvnrnd(x_t_mean, h^2/M * eye(p));
end
fprintf('Sampling done (finished in %.1f seconds). Acceptance rate=%f\n', toc, 1-count_rej/M/N);
end
