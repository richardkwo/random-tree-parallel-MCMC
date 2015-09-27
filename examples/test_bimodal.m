addpath ../src/PART ../src/PART/cut ../src/alternatives ../src/utils

%% Combine densities with two modes
% The data is generated using the following generative model: 
%
% $$x \sim \prod_{i = 1}^m p_i(x)$$
%
% where m is the number of subsets, n is the number of observations and
% each $p_i(\cdot)$ is a mixture of two Gaussians,
%
% $$ p_i(x) = w_{i1} N(u_{i1}, s_{i1}^2) + w_{i2} N(u_{i2}, s_{i2}^2)$$

p = 1; % dim
m = 10; % number of subsets
n = 10000; % number of total data
u = [0 + normrnd(-5,0.5,m,1), 5 + normrnd(0,0.5,m,1)]; % means
s = [1 + abs(normrnd(0,0.1,m,1)), 4 + abs(normrnd(0,0.1,m,1))]; %standard deviations
w = zeros(m,2); %weights
w(:,1) = 0.8/3;
w(:,2) = 1 - w(:,1);

% Here we show the true density of the data
x = -20:0.001:20;
y = zeros(m, length(x));
for i = 1:length(x)
y(:,i) = w(:,1).*normpdf(x(i), u(:,1), s(:,1)) + w(:,2).*normpdf(x(i), u(:,2), s(:,2));
end

y1 = prod(y);
y1 = y1/sum(y1)/0.001;

figure;
ymax = max(y1) + 0.05;
axis_f = [-10 10 0 ymax];
plot(x,y1,'-r');
axis(axis_f)
hold on;
title('The true density for the multiplied distribution')

N = 5000; %Number of samples drawn from each subset.
posterior_N = N;
sub_chains = cell(1,m);
sub_density = cell(1,m);
for i = 1:m
    sub_chains{i} = [normrnd(u(i,1), s(i,1), floor(w(i,1)*N), 1);...
        normrnd(u(i,2), s(i,2), N - floor(w(i,1)*N), 1)];
    sub_chains{i} = sub_chains{i}(randperm(size(sub_chains{i},1)),:);
    % subset density
    sub_density{i} = w(i,1) * normpdf(x, u(i,1), s(i,1)) + w(i,2) * normpdf(x, u(i,2), s(i,2));
    plot(x, sub_density{i}, '--');
end


%% Combine via PART: one-stage aggregation
options = part_options('min_cut_length', 0.001, 'min_fraction_block', 0.01, 'ntree', 8, 'verbose', 2);
combined_posterior_kd_onestage = aggregate_PART_onestage(sub_chains, options);
combined_posterior_ml_onestage = aggregate_PART_onestage(sub_chains, options);
options.local_gaussian_smoothing = false;
combined_posterior_kd_onestage_NoVar = aggregate_PART_onestage(sub_chains, options);
combined_posterior_ml_onestage_NoVar = aggregate_PART_onestage(sub_chains, options);


figure;
plot_marginal_compare({combined_posterior_kd_onestage, ...
    combined_posterior_ml_onestage, ...
    combined_posterior_kd_onestage_NoVar, ...
    combined_posterior_ml_onestage_NoVar}, ...
    { 'KD-onestage', 'ML-onestage','KD-onestage-NoSmoothing','ML-onestage-NoSmoothing'});
subplot(1,1,1);
hold on;
plot(x,y1,'--','DisplayName', 'True', 'LineWidth',2);
legend('-DynamicLegend', 'Location', 'best'); 
title('PART: one-stage aggregation'); hold off;


%% Combine via PART: pairwise aggregation
options = part_options('min_cut_length', 0.001, 'min_fraction_block', 0.01, 'ntree', 8, 'verbose', 2);
combined_posterior_kd_pairwise = aggregate_PART_pairwise(sub_chains, options);
combined_posterior_ml_pairwise = aggregate_PART_pairwise(sub_chains, options);

figure;
plot_marginal_compare({combined_posterior_kd_pairwise, ...
    combined_posterior_ml_pairwise}, ...
    { 'KD-pairwise','ML-pairwise'});
subplot(1,1,1);
hold on;
plot(x,y1,'--','DisplayName', 'True', 'LineWidth',2);
legend('-DynamicLegend', 'Location', 'best'); 
title('PART: pairwise (multi-stage) aggregation'); hold off;


%% Other methods for comparison

% simple averaging
combined_posterior_averaging = aggregate_average(sub_chains);

% weighted averaging 
combined_posterior_weighted_averaging = aggregate_weighted_average(sub_chains);

% UAI - parametric
combined_posterior_parametric = aggregate_uai_parametric(sub_chains);

% UAI - nonparametric
combined_posterior_nonparametric = aggregate_uai_nonparametric(sub_chains, 1e4);

% UAI - semiparametric
combined_posterior_semiparametric = aggregate_uai_semiparametric(sub_chains, 1e4);

figure;
plot_marginal_compare({combined_posterior_kd_pairwise, combined_posterior_ml_pairwise, ...
    combined_posterior_averaging, ...
    combined_posterior_weighted_averaging, ...
    combined_posterior_parametric, ...
    combined_posterior_nonparametric, ...
    combined_posterior_semiparametric}, ...
    {'PART-KD', 'PART-ML', ...
    'average', 'weighted average', ...
    'Neiswanger - parametric', 'Neiswanger - nonparametric', 'Neiswanger - semiparametric'});
subplot(1,1,1);
hold on;
plot(x,y1,'--','DisplayName', 'True', 'LineWidth',2);
legend('-DynamicLegend', 'Location', 'best'); 
title('Comparison of various algorithms'); hold off;
