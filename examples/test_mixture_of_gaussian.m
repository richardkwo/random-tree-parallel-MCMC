addpath ../src/PART ../src/PART/cut ../src/alternatives ../src/utils

%% generate mixture of mv-normal sample, 2-D
K = 2; % number of components
M = 4; % split
n = 400;
Sigmas = cell(1,K);
Mus = zeros(K, 2);

% draw means
mean_mu = [5, 5];
sigma_mu = 10 * eye(2);
Mus = mvnrnd(repmat(mean_mu, K, 1), sigma_mu);

% use the same cov for now
for k=1:K
    Sigmas{k} = eye(2);
end

% mixture weights
mixture_weights = dirichrnd(1, 1 * ones(1,K));

% draw samples
all_data = zeros(n, 2);
for i=1:n
    z =randsample(K, 1, true, mixture_weights);
    all_data(i,:) = mvnrnd(Mus(z,:), Sigmas{z});
end

% split into subsets of data
subset_data = cell(1,M);
assign_vec = repmat(1:M, 1, n/M);
for c=1:M
    subset_data{c} = all_data(assign_vec==c, :);
end

rr = floor(sqrt(M+1)); ss = ceil(sqrt(M+1));
for c=0:M
	if c==0
		figure;
		subplot(rr,ss,1);
		plot(all_data(:,1), all_data(:,2), 'k.');
		hold on;
		plot(Mus(:,1)', Mus(:,2)', 'ro');
		title('full data');
	else
		subplot(rr,ss,c+1);
		plot(subset_data{c}(:,1), subset_data{c}(:,2), 'k.');
		hold on;
		plot(Mus(:,1)', Mus(:,2)', 'ro');
		title(sprintf('subset data %d', c));
	end
end

%% mcmc

N = 5000;
fprintf('Sampling full chain...\n');
full_chain = sample_mog(all_data, Sigmas, Mus, N);
sub_chain = cell(1, M);
for c=1:M
    fprintf('Sampling subset %d...\n', c);
    
    sub_chain{c} = sample_mog(subset_data{c}, Sigmas, Mus, N);
end

rr = floor(sqrt(M+1)); ss = ceil(sqrt(M+1));
for c=0:M
	if c==0
		figure;
		subplot(rr,ss,1);
		plot(full_chain(:,1), full_chain(:,2), 'k.');
		hold on;
		plot(Mus(:,1)', Mus(:,2)', 'ro');
		title('full chain');
	else
		subplot(rr,ss,c+1);
		plot(sub_chain{c}(:,1), sub_chain{c}(:,2), 'k.');
		hold on;
		plot(Mus(:,1)', Mus(:,2)', 'ro');
		title(sprintf('subchain %d', c));
	end
end



%% averaging, weighted averaging & parametric
combined_posterior_averaging = aggregate_average(sub_chain);
combined_posterior_weighted_averaging = aggregate_weighted_average(sub_chain);
combined_posterior_parametric = aggregate_uai_parametric(sub_chain);

figure;
plot(combined_posterior_averaging(:,1), combined_posterior_averaging(:,2), 'rx'); hold on;
plot(combined_posterior_weighted_averaging(:,1), combined_posterior_weighted_averaging(:,2), 'bo');
plot(combined_posterior_parametric(:,1), combined_posterior_parametric(:,2), 'gs');
plot(full_chain(:,1), full_chain(:,2), 'k.'); 
legend('average', 'weighted average', 'parametric', 'True');

%% KD/ML aggregations

options = part_options('min_cut_length', 0.01, 'min_fraction_block', 0.01, 'resample_N', 5000, 'ntree', 16);
combined_posterior_kd_pairwise = aggregate_PART_pairwise(sub_chain, options);
options.cut_type = 'ml';
[combined_posterior_ml_pairwise, sampler, ~] = aggregate_PART_pairwise(sub_chain, options);

figure;
plot(combined_posterior_kd_pairwise(:,1), combined_posterior_kd_pairwise(:,2), 'rx'); hold on;
plot(combined_posterior_ml_pairwise(:,1), combined_posterior_ml_pairwise(:,2), 'bo');
plot(full_chain(:,1), full_chain(:,2), 'k.'); 
legend('PART-KD', 'PART-ML', 'True', 'Location', 'best');

figure;
for i=1:16
subplot(4,4,i); axis off;
plot_tree_blocks(sampler{i}, []);
title(['Tree ', num2str(i)]);
end

%% other combinations
N = 1e5; ii = randsample(N, 1000);
% UAI - nonparametric
combined_posterior_nonparametric = aggregate_uai_nonparametric(sub_chain, N);
% UAI - semiparametric
combined_posterior_semiparametric = aggregate_uai_semiparametric(sub_chain, N);

figure;
plot(combined_posterior_nonparametric(ii,1), combined_posterior_nonparametric(ii,2), 'rx'); hold on;
plot(combined_posterior_semiparametric(ii,1), combined_posterior_semiparametric(ii,2), 'bo');
plot(full_chain(:,1), full_chain(:,2), 'k.'); 
legend('nonparametric', 'semiparametric', 'True');

%% summary
performance_table({full_chain, ...
    combined_posterior_kd_pairwise, ...
    combined_posterior_ml_pairwise, ...
    combined_posterior_averaging, ...
    combined_posterior_weighted_averaging, ...
    combined_posterior_parametric, ...
    combined_posterior_nonparametric, ...
    combined_posterior_semiparametric}, ...
    {'true', ...
    'KD', 'ML', ...
    'average', 'weighted average', ...
    'Neiswanger - parametric', 'Neiswanger - nonparametric', 'Neiswanger - semiparametric'}, ...
    full_chain)
fprintf('NOTE: this table may not be a good evaluation since the posterior is non-Gaussian.\n');
