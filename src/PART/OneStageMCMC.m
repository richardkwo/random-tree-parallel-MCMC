function [trees, sampler, sampler_prob, RawMCMC] = OneStageMCMC(MCdraws, option)
% OneStageMCMC is a function that takes an input of raw subset posteriors
% and package it into a RawMCMC structure. It will also call buildForest
% automatically to produce the combined densities in an one-stage manner
% 
% Calls:
% [trees, sampler, sampler_prob] = OneStageMCMC(MCdraws, option)
%   
% Arguments:
% 
% MCdraws: a 1 x m cell containing subset posteriors from m subsets
% option: option configured by part_options(...). Use `option = part_options()` for default settings. 
% 
% Outputs:
% 
% trees: A list containing multiple partition trees, each corresponding
%        to a density estimation of the combined posterior. Used by
%        treeDensity to evaluate density at a given point
% sampler: A list containing multiple flattened partition trees, used by
%          treeSampling to resample from the forest
% sampler_prob: A list containing the probabilities for each flattened
%               partition tree, used be called by treeSampling
% 
% See also:
% buildForest, MultiStageMCMC, NormalizeTree, treeSampling, treeDensity
% 
    
    % setting up the default value for parameters.

    if ~isfield(option,'area'), area = [];else area = option.area;end   
    if ~isfield(option, 'rule')
        switch option.cut_type
            case 'kd'
                option.rule = @kdCut;
            case 'ml'
                option.rule = @consensusMLCut;
        end
    end 
    m = length(MCdraws); %number of subsets
    d = size(MCdraws{1}, 2); %number of dimensions
    N = zeros(1, m); %Store number of posterior samples for different subsets
    
    %Initialization
    list = []; %To concatenate MCMC samples
    mark = []; %To index different points
    
    for i = 1:m
        N(i) = size(MCdraws{i},1); %Obtain the subset posterior size
        list = [list;MCdraws{i}]; %Concatenate the posterior draws
        mark = [mark;[1:N(i);ones(1,N(i))*i]']; %augmenting the index list
    end

    if ~isfield(option, 'min_number_points_block')
        RawMCMC.min_number_points_block = ceil(max(N) * option.min_fraction_block);
    else 
        RawMCMC.min_number_points_block = option.min_number_points_block; 
    end

    
    %Packaging all information into RawMCMC
    RawMCMC.list = list;
    RawMCMC.area = area;
    RawMCMC.N = N;
    RawMCMC.para = 1:d;
    RawMCMC.mark = mark;
    RawMCMC.total_set = m;
    RawMCMC.ntree = option.ntree;
    RawMCMC.rule = option.rule;
    RawMCMC.parallel = option.parallel;
    RawMCMC.verbose = option.verbose;
    RawMCMC.min_cut_length = option.min_cut_length;
    
    %Call buildForest to produce the combined density
    fprintf('Building tree ensemble: minimum %d points per block, minimum side length = %f\n', RawMCMC.min_number_points_block, RawMCMC.min_cut_length);
    [trees, sampler, sampler_prob] = buildForest(RawMCMC, option.verbose);
    block_numbers = zeros(1, option.ntree);
    for t=1:option.ntree
        block_numbers(t) = length(sampler{t});
    end
    fprintf('\nRandom tree ensemble constructed -- %d trees and %f leafs per tree\n', option.ntree, mean(block_numbers));
end