function re = relative_error(samples, full_chain, theta)
% re = relative_error(samples, full_chain, theta)
% ||posterior_theta - theta||_2 / ||posterior_theta (full chain) - theta||_2
theta = theta(:)';
assert(size(samples,2)==length(theta) && size(full_chain,2)==length(theta));
error_full_chain = mean((full_chain - repmat(theta, size(full_chain,1), 1)).^2, 1);
error_samples = mean((samples - repmat(theta, size(samples, 1), 1)).^2, 1);
re = sqrt(sum(error_samples) / sum(error_full_chain));
end

