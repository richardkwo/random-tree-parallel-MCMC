function rmse = rmse_posterior_mean(samples_left, samples_right)
% rmse of the posterior mean of two samples
assert(size(samples_left,2)==size(samples_right,2),'dim not match');
mean_left = mean(samples_left,1);
mean_right = mean(samples_right,1);
rmse = sqrt(mean((mean_left - mean_right).^2));
end

