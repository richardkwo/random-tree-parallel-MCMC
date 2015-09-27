function [newForest, sampler, sampler_prob] = NormalizeForest(Forest, MCdraws)
    %NormalizeForest is a function to process a forest of unprocessed
    %random partition trees. The situation that we need this function is
    %that when a tree stucture is given, and we insert subset posteriors into
    %that structure, then we need process the forest to obtain probability,
    %density, node.mean and node.var and output the corresponding
    %partition tree and sampler. NormalizeForest calls NormalizeTree to
    %normalize each of its component and produce the corresponding sampler.
    %
    %Calls:
    %[newForest, sampler, sampler_prob] = NormalizeForest(Forest, MCdraws)
    %
    %Arguments:
    %
    %Forest: A forest of unprocessed partition tree. Typically resulted
    %   from InsertNode.m
    %MCdraws: A list of raw MCMC samples from subsets
    %
    %Output:
    %
    %newForest: Normalized Forest that can be used by treeDensity
    %sampler: Flattened normalized forest that can be used by treeSampling
    %sampler_prob: corresponding log(probability) of sampler
    %
    %See also:
    %NormalizeTree
    %
    
    ntree = length(Forest); %Number of trees
    
    %Initialization
    newForest = cell(1,ntree); 
    sampler = cell(1,ntree);
    sampler_prob = cell(1,ntree); 
    
    %Setting up waitbar
    h = waitbar(0, 'Starting normalizing forest...');
    
    %Calling NormalizeTree to process
    for i = 1:ntree
        waitbar((i-1)/ntree, h, ['Evaluating the ', num2str(i), 'th tree...']);
        [a, b, p] = NormalizeTree(Forest{i}, MCdraws);
        newForest{i} = a;
        sampler{i} = b;
        sampler_prob{i} = p;
    end
    close(h);
end