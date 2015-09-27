function [aggregated_samples, sampler_onestage, prob_onestage] = aggregate_PART_onestage(sub_chains, options)
% Aggregating sub-chain samples with PART by one-stage combination
% 
% Call:
% [aggregated_samples, sampler_onestage, prob_onestage] = aggregate_PART_onestage(sub_chains, options)
% 
% Arguments: 
% sub_chains: a 1 x m cell of m sub-chains. 
%             sub_chains{i} should be an N_i x p matrix for N_i number of p-dimensional MCMC samples.  
% options:    options configured by part_options(...). Use options = part_options() for default settings.
% 
% Output: 
% aggregated_samples: aggregated posterior samples. Number of samples set by options.resample_N.
% sampler_onestage: flattened partition trees, used by is used by treeSampling(...) and teeDensity(...).
% prob_onestage: normalized leaf probabilities corresponding to sampler_onestage, used by is used by treeSampling(...) and teeDensity(...).

N = max([size(sub_chains{1}, 1), options.resample_N]);
M = length(sub_chains);

options.UseData = false;
options.list = [];
options.mark = [];

if strcmp(options.cut_type, 'ml')
    options.rule = @consensusMLCut;
elseif strcmp(options.cut_type, 'kd')
    options.rule = @kdCut;
end

[~, sampler_onestage, prob_onestage] = OneStageMCMC(sub_chains, options);
aggregated_samples = treeSampling(sampler_onestage, prob_onestage, N, options);

end

