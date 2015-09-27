function T = performance_table(aggregated_posteriors, labels, full_posterior, theta)
L = length(aggregated_posteriors);
assert(L==length(labels));
metric_names = {'rmse_cov', 'rmse_var', 'rmse_mean', 'KL_true_vs_est', 'KL_est_vs_true', 'relative_l2_error'};
metrics = zeros(L, length(metric_names));
for l=1:L
    for i=1:length(metric_names)
        switch metric_names{i}
            case 'rmse_cov'
                metrics(l, i) = rmse_posterior_cov(aggregated_posteriors{l}, full_posterior);
            case 'rmse_var'
                metrics(l, i) = rmse_posterior_cov(aggregated_posteriors{l}, full_posterior, 'diag');
            case 'rmse_mean'
                metrics(l, i) = rmse_posterior_mean(aggregated_posteriors{l}, full_posterior);
            case 'KL_true_vs_est'
                metrics(l, i) = approximate_KL(full_posterior, aggregated_posteriors{l});
            case 'KL_est_vs_true'
                metrics(l, i) = approximate_KL(aggregated_posteriors{l}, full_posterior);
            case 'relative_l2_error'
                metrics(l, i) = relative_error(aggregated_posteriors{l}, full_posterior, theta);
        end
    end
end
T = array2table(metrics, 'VariableNames', metric_names, 'RowNames', labels);
end

