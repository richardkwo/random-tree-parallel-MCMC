function p = dirichrnd(n, alpha)
    % a density of prod(x_i^(alpha_i - 1))
    m = length(alpha);
    internal = gamrnd(repmat(alpha, n, 1),1);
    p = internal./repmat(sum(internal, 2),1,m);
end