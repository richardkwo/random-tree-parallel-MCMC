function aggregated_samples = aggregate_uai_parametric(sub_chain, varargin)
% "parametric" aggregation from the UAI paper
% multiplicated gaussian with laplacian approximation to each subchain
if ~isempty(varargin)
    n = varargin{1};
else
    n = size(sub_chain{1},1);
end
M = length(sub_chain);
[~,p] = size(sub_chain{1});
Prec = cell(1,M);
agg_Sig = zeros(p);
agg_Mu = zeros(1,p);
for c=1:M
    Prec{c} = inv(cov(sub_chain{c}));
    agg_Sig = agg_Sig + Prec{c};
    agg_Mu  = agg_Mu +  mean(sub_chain{c},1) * Prec{c};
end
agg_Sig = inv(agg_Sig);
agg_Mu = agg_Mu * agg_Sig;

aggregated_samples = mvnrnd(repmat(agg_Mu,n,1), agg_Sig);
end