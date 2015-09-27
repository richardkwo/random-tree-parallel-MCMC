function rmse = rmse_posterior_cov(samples_left, samples_right, varargin)
% compute the rmse of posterior cov between two samples
assert(size(samples_left,2)==size(samples_right,2),'dim not match');
cov_left = cov(samples_left);
cov_right = cov(samples_right);
errors = (cov_left - cov_right).^2;
if ~isempty(varargin) && strcmp(varargin{1},'diag')
    errors = diag(errors);
end
rmse = sqrt(mean(errors(:)));
end