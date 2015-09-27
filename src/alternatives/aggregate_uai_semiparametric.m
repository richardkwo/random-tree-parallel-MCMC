function aggregated_samples = aggregate_uai_semiparametric(sub_chain, varargin)
% Drawing N samples from aggregated posterior with Neiswager et al's 
% semiparametric method

if ~isempty(varargin)
    N = varargin{1};
else
    N = size(sub_chain{1},1);
end

M = length(sub_chain);
[~,p] = size(sub_chain{1});
subchain_sizes = zeros(1,M);
for c=1:M
    subchain_sizes(c) = size(sub_chain{c},1);
end

subchain_Sigma = cell(1,M);
subchain_Mu = cell(1,M);
Prec = cell(1,M);
sum_of_Prec = zeros(p);
agg_Mu = zeros(1,p);
for c=1:M
    subchain_Sigma{c} = cov(sub_chain{c});
    subchain_Mu{c} = mean(sub_chain{c},1);
    Prec{c} = inv(subchain_Sigma{c});
    sum_of_Prec = sum_of_Prec + Prec{c};
    agg_Mu  = agg_Mu + subchain_Mu{c} * Prec{c};
end
agg_Sigma = inv(sum_of_Prec);
agg_Mu = agg_Mu / sum_of_Prec;

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
whole_normal = log_pdf_normal_whole(x_t_mean, agg_Mu, agg_Sigma, 1, M, p);
subchain_normals = log_pdf_subchain_normals(x_t, subchain_Mu, subchain_Sigma, M);

fprintf('Getting %d samples with semi-parametric method...\n', N);
tic
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
        whole_normal_new = log_pdf_normal_whole(x_t_new_mean, agg_Mu, agg_Sigma, h, M, p);
        subchain_normals_new = log_pdf_subchain_normals(x_t_new, subchain_Mu, subchain_Sigma, M);
        u = rand;
        if u < exp(-(ss_new - ss_t)/(2*h^2) + ...
                whole_normal_new - whole_normal + subchain_normals - subchain_normals_new)
            % accept the move
            index_t = index_new;
            x_t = x_t_new;
            x_t_mean = x_t_new_mean;
            ss_t = ss_new;
            whole_normal = whole_normal_new;
            subchain_normals = subchain_normals_new;
        else
            count_rej = count_rej + 1;
        end
    end
    % sample with the index
    aggregated_samples(i, :) = mvnrnd(x_t_mean, h^2/M * eye(p));
end
fprintf('Sampling done (finished in %.1f seconds). Acceptance rate=%f\n', toc, 1-count_rej/M/N);

end

function l = log_pdf_normal_whole(x_t_mean, agg_Mu, agg_Sig, h, M, p)
l = logmvnpdf(x_t_mean, agg_Mu, agg_Sig + h/M*eye(p));
end

function l = log_pdf_subchain_normals(x_t, subchain_Mu, subchain_Sigma, M)
l = 0;
for c=1:M
    l = l + logmvnpdf(x_t(c,:), subchain_Mu{c}, subchain_Sigma{c});
end
end
