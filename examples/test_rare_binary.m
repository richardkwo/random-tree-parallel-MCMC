addpath ../src/PART ../src/PART/cut ../src/alternatives ../src/utils

%% Binomial data with rare event
% The data is generated using the following generative model: 
%
% $$x \sim bin(p), \quad p = 2m/n$$
%
% where m is the number of subsets, n is the number of observations.

p = 1; % 1 dim
cat = 2; %binomial distribution
M = 15; %number of subsets
n = 10000; %number of observations
theta = zeros(1, cat); %probability
theta(1) = 2*M/n; 
theta(2) = 1 - theta(1);

% beta prior on p
pseudocounts = [2 2];

N = 10000; % number drawn for posterior
X = mnrnd(1, theta, n);

% We then random partition the data set into m subsets
sub_X = cell(1,M);
for m = 1:M
    sub_X{m} = X(m:M:end, :);
end
subset_totals = cellfun(@(x) sum(x(:,1)), sub_X)

% We draw N posterior samples for full data
full_chain = dirichrnd(N, sum(X, 1) + pseudocounts);
full_chain = full_chain(:,1);

% We draw N posterir samples for each subset posterior
sub_chain = cell(1,M);
for m = 1:M
    sub_chain{m} = dirichrnd(N, sum(sub_X{m}, 1) + (pseudocounts/M + (M-1)/M));
    sub_chain{m} = sub_chain{m}(:,1);
end


do_posterior_comparison = true;

if do_posterior_comparison
    fprintf('Plotting full chain and sub-chains...\n');
    figure;           
    hold on;
    [f,xi] = ksdensity(full_chain);
    plot(xi,f,'DisplayName','full');
    for l=1:M
        [f,xi] = ksdensity(sub_chain{l});
        plot(xi,f,'DisplayName',['subset ',num2str(l)]);
    end        
%     legend('-DynamicLegend', 'Location', 'best'); 
%     title('theta');
    hold off;

end

%% parametric
combined_posterior_averaging = aggregate_average(sub_chain);
combined_posterior_weighted_averaging = aggregate_weighted_average(sub_chain);
combined_posterior_parametric = aggregate_uai_parametric(sub_chain);
figure;
plot_marginal_compare({full_chain, combined_posterior_parametric, combined_posterior_averaging, combined_posterior_weighted_averaging}, ...
    {'true', 'parametric','average','weighted'});

%% KD/ML aggregations

options = part_options('min_cut_length', 0.01, 'resample_N', 10000);
combined_posterior_kd_pairwise = aggregate_PART_pairwise(sub_chain, options);
options.cut_type = 'ml';
combined_posterior_ml_pairwise = aggregate_PART_pairwise(sub_chain, options);

figure;
plot_marginal_compare({full_chain, ...
    combined_posterior_kd_pairwise, combined_posterior_ml_pairwise}, ...
    {'true', 'kd-pairwise', 'ml-pairwise'});
xlim([-2 14]*1e-3);


%% other combinations
% UAI - nonparametric
combined_posterior_nonparametric = aggregate_uai_nonparametric(sub_chain, 1e4);

% UAI - semiparametric
combined_posterior_semiparametric = aggregate_uai_semiparametric(sub_chain, 1e4);

%% plotting
% marginal density plot
figure;
plot_marginal_compare({full_chain, sub_chain{1} ,...
    combined_posterior_kd_pairwise, ...
    combined_posterior_ml_pairwise, ...
    combined_posterior_averaging, ...
    combined_posterior_weighted_averaging, ...
    combined_posterior_parametric, ...
    combined_posterior_nonparametric, ...
    combined_posterior_semiparametric}, ...
    {'true', 'subset 1',...
    'KD', 'ML', ...
    'average', 'weighted average', ...
    'Neiswanger - parametric', 'Neiswanger - nonparametric', 'Neiswanger - semiparametric'});
xlim([0 15]*1e-3);

% comparison of accuracy
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
    full_chain, theta(1))