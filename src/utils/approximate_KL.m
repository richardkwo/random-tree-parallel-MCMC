function kl = approximate_KL(samples_left, samples_right)
% compute the apprximate KL (p_left || p_right)
% The KL is computed from Laplacian approximations fitted to both
% distributions
assert(size(samples_left,2)==size(samples_right,2),'dim not match');
p = size(samples_left,2);
mu_left = mean(samples_left,1);
mu_right = mean(samples_right, 1);
Sigma_left = cov(samples_left);
Sigma_right = cov(samples_right);
inv_Sigma_right = inv(Sigma_right);
d_mu = mu_right - mu_left;
kl = 1/2 * (trace(inv_Sigma_right * Sigma_left) + ...
    d_mu * inv_Sigma_right * d_mu' - p + log(det(Sigma_right)) - log(det(Sigma_left)));
end

