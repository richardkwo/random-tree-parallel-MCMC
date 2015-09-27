function ll = empiricalLogLikelihood(x, l, r, cuts)
% ll = empiricalLogLikelihood(x, l, r, cuts)
% evaluates the empirlcal log likelihood of x under cuts
%
% x: data
% l: left boundary
% r: right boundary
% cuts: cutting points
% 
% ll: empirical log likelihood
assert(all(x>=l) && all(x<=r), 'not all x in the range [l,r]');
cuts = [l sort(cuts) r];
N = length(x);
ll = 0;
for k=1:(length(cuts)-1)
    lk = cuts(k);
    rk = cuts(k+1);
    nk = sum(lk<=x & x<rk);
    ll = ll + nk * log (nk / N / (rk-lk)); % ll = sum_k n_k * log ((n_k/N)/|A_k|)
end
end