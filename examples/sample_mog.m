function mu = sample_mog(data, Sigmas, Mus, N)
n = size(data, 1);
K = length(Sigmas);
p = size(Mus, 2);
mu = zeros(N, p);
Z = zeros(n, K);
W = zeros(n, K);
for i=1:N
	% sample assignments
	for k=1:K
		W(:,k) = mvnpdf(data, Mus(k,:), Sigmas{k});
	end
	for j=1:n
		Z(j, :) = 0;
		Z(j, randsample(K, 1, true, W(j,:))) = 1;
	end
	% sample mean
	for k=1:K
		Mus(k,:) = mvnrnd(mean(data(Z(:,k)==1, :)), Sigmas{k}/sum(Z(:,k)));
	end
	kk = randsample(K, 1, true, sum(Z));
	mu(i,:) = Mus(kk, :);
end
mu = mu(randperm(N), :);
end