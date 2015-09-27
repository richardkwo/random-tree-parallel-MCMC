function [trees, sampler, sampler_prob] = buildForest(RawMCMC, verbose)
    %buildForest is a function that takes input of a list of one or multiple subset posteriors
    %to produce the combined posterior density using an ENSEMBLE of random partition trees. 
    %The function calls buildTree multiple times to build a number (specified by user) 
    %of iid partition trees. This function can be used alone or
    %called by CombineMCMC and MultiMCMC for combining with different
    %strategies
    %
    %Calls:
    %[trees, sampler, sampler_prob] = buildForest(RawMCMC, verbose)
    %
    %Arguments:
    %
    %RawMCMC is a data structure containing:
    %   list: a list of MCMC samples from different subsets
    %   area: the total area for cutting, the
    %   N: number of posterior samples of each subset MCMC samples
    %   para: the columns/dimensions that represent parameters of interest
    %   mark: Index of data point. It would cost too much to store data directly
    %      in "obj" and "p", we instead store their index within each node
    %      in a form of (i, j) where j is the index for the subset and i is
    %      the index of the data point in that subset. For example, (1, 2) 
    %      means the first point from the second subset
    %   total_set: total number of subsets
    %   ntree: the number of random partition trees to grow
    %   rule: A function handle that determines the optimal cut at each gate
    %   parallel: if True, building trees in a forest in parallel
    %   max_n_unique: FORCE it to be a leaf node if has < max_n_unique number of unique samples in a region
    %   min_number_points_block: stopping criteria: continue cutting if the number of points in a block > min_number_points_block
    %   min_cut_length: either a scalar of a vector of the same dim as parameters. The lengths of a block >= min_cut_length
    %
    %verbose: indicateing whether a waitbar should be shown
    %
    %Outputs:
    %
    %trees: A list containing multiple partiiton trees, each corresponding
    %       to a density estimation of the combined posterior. Used by
    %       treeDensity to evaluate density at a given point
    %sampler: A list containing multiple flattened partition trees, used by
    %         treeSampling to resample from the forest: 
    %         {(leafs/blocks of tree_1), (leafs/blocks of tree_2), ...}
    %sampler_prob: A list containing the probabilities for each flattened
    %              partition tree, used be called by treeSampling:
    %              {(probs for blocks of tree_1), (probs for blocks of tree_2), ...}
    %
    %See also:
    %buildTree, OneStageMCMC, MultiStageMCMC, NormalizeTree, treeSampling,
    %treeDensity
    %
    
    %By default, we use a waitbar
    if nargin == 1
        verbose = 0;
    end
    
    %Obtain the list of subset posteriors
    if ~isfield(RawMCMC, 'list'), error('no data loaded.');else list = RawMCMC.list; end
    [n, d] = size(list);
    
    %Check the different component of RawMCMC and specify the default value
    if ~isfield(RawMCMC, 'area'), area = []; else area = RawMCMC.area; end
    if ~isfield(RawMCMC, 'N'), error('I dont know the number of posterior samples on each subset.'), else N = RawMCMC.N; end
    if ~isfield(RawMCMC, 'para'), para = 1:d; else para = RawMCMC.para; end
    if ~isfield(RawMCMC, 'mark'), mark = [1:n;ones(1,n)]'; else mark = RawMCMC.mark; end
    if ~isfield(RawMCMC, 'total_set'), total_set = 1; else total_set = RawMCMC.total_set; end
    if ~isfield(RawMCMC, 'ntree'), ntree = 20; else ntree = RawMCMC.ntree; end
    if ~isfield(RawMCMC, 'rule'), rule = @kdCut; else rule = RawMCMC.rule; end
    if ~isfield(RawMCMC, 'parallel'), do_parallel=false; else do_parallel=RawMCMC.parallel; end
    if ~isfield(RawMCMC, 'max_n_unique'), max_n_unique=10; else max_n_unique=RawMCMC.max_n_unique; end
    if ~isfield(RawMCMC, 'min_number_points_block'), min_number_points_block=20; else min_number_points_block=RawMCMC.min_number_points_block; end
    if ~isfield(RawMCMC, 'min_cut_length'), min_cut_length=0.1; else min_cut_length=RawMCMC.min_cut_length; end

    %Initialization. Cell lists to restore information for different random
    %partition trees
    trees = cell(1,ntree);
    sampler = cell(1,ntree); 
    sampler_prob = cell(1,ntree);
    
    %Set up waitbar
    if verbose>1
        h = waitbar(0,'The final stage...');
    end
    
    parfor_progress(ntree);
    %Start building forest...
    if do_parallel && ntree>1
        fprintf('\nBuilding a forest of %d trees in parallel...\n', ntree);
        parfor i = 1:ntree
            %Call buildTree to build random partition tree
            [tree, nodes, log_probs] = buildTree(list, area, N, para, mark, total_set, ...
                min_number_points_block, min_cut_length, 0, rule, max_n_unique);
            % log of normalizer
            log_normalizer = 0;
            %When there are multiple subsets, we need an extra normalizing step
            if total_set > 1
                %normalizing densities using the total probability
                % first substracting the max to prevent enormous exp(..)
                max_log_prob = max(log_probs);
                log_normalizer = max_log_prob + log(sum(exp(log_probs - max_log_prob)));
                tree = NormalizeTree(tree, [], log_normalizer); 
            end

            %save the results
            trees{i} = tree;
            sampler{i} = nodes;
            sampler_prob{i} = log_probs - log_normalizer; % normalized probs in log scale
            
            assert(abs(sum(exp(sampler_prob{i}))-1)<1e-10) % FIXME: comment it out
            
            % text update (only this works for parallel)
            parfor_progress;
        end
        parfor_progress(0);
    else
        fprintf('\nBuilding a forest of %d trees...\n', ntree);
        for i = 1:ntree
            fprintf('Building %d-th out of %d trees...\n', i, ntree);
            %Call buildTree to build random partition tree
            [tree, nodes, log_probs] = buildTree(list, area, N, para, mark, total_set, ...
                min_number_points_block, min_cut_length, 0, rule, max_n_unique);
            % log of normalizer
            log_normalizer = 0;
            %When there are multiple subsets, we need an extra normalizing step
            if total_set > 1
                %normalizing densities using the total probability
                % first substracting the max to prevent enormous exp(..)
                max_log_prob = max(log_probs);
                log_normalizer = max_log_prob + log(sum(exp(log_probs - max_log_prob)));
                tree = NormalizeTree(tree, [], log_normalizer); 
            end

            %save the results
            trees{i} = tree;
            sampler{i} = nodes;
            sampler_prob{i} = log_probs - log_normalizer; % normalized probs in log scale
            
            assert(abs(sum(exp(sampler_prob{i}))-1)<1e-10) % FIXME: comment it out

            % update waitbar
            if verbose>1
                waitbar(i/ntree,h);
            end
        end
    end
    %close waitbar
    if verbose>1
        close(h);
    end
end