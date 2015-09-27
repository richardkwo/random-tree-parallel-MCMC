function y = treeSampling(sampler, sampler_prob, n, option)
    %treeSampling use the output from buildForest, OneStageMCMC and
    %MultiStageMCMC to resample points from the combined density.
    %
    %Call:
    %y = treeSampling(sampler, sampler_prob)
    %y = treeSampling(sampler, sampler_prob, n)
    %y = treeSampling(sampler, sampler_prob, n, option)
    %
    %Arguments:
    %
    %sampler: the output from buildForest, OneStageMCMC and MultiStageMCMC,
    %   which is a forest of flattened random partition tree
    %   sampler_prob: the corresponding probiblity for sampler, also outputed
    %   by the three functions aforementioned
    %n: number of resamples
    %option: optional choices including:
    %   resample_data_only: an indicator whether the resampling within each stage
    %       should use the original data points, i.e., when resampling at a
    %       given leaf, we only resample original data points that are within
    %       that leaf.
    %   local_gaussian_smoothing: an indicator whether the resampling should use an normal
    %       density on each leaf, i.e., when resampling at a given leaf,
    %       we use the original data points located in the leaf to fit a
    %       normal density (the variance will be shrinked by a factor of
    %       1/total_subset), and resample using this density
    %   list: the list containing original subset posteriors concatenated
    %       in one matrix (See buildForest)
    %   mark: index for original subset posteriors (See buildForest)
    %
    %Outputs:
    %
    %y: the resampled posterior samples from the combined density
    %
    %See also:
    %buildForest, OneStageMCMC, MultiStageMCMC, treeDensity
    %
    
    %The default number of resampling
    if nargin == 2
        n = 1;
    
    %The default options
    elseif nargin == 3
        option.list = [];
        option.mark = [];
        option.resample_data_only = false;
        option.local_gaussian_smoothing = false;
    elseif nargin == 4
        %Check if the option is consistent
        if ~isempty(option.list) && isempty(option.mark)
            error('mark must be provided if list exists');
        elseif option.resample_data_only && isempty(option.list)
            error('Need to provide the data in order to resample_data_only');
        end
    end
    
    m = length(sampler); %number of trees
    d = size(sampler{1}{1}.area,1); %data dimension
    y = zeros(n, d); %restore the resamples
    
    %h = waitbar(0,['Sampling ',num2str(n), ' %samplers from the distribution...']);
    
    if option.resample_data_only
        %if resampling using the original subset posteriors, we need to
        %obtain the sizes of posterior samples on different subsets
        N = zeros(1,length(unique(option.mark(:,2))));
        for set = 1:length(N)
            N(set) = sum(option.mark(:,2)< set);
        end
        %N is accumulated sizes i.e., N(1) = 0, N(2) = |subset1|, N(3) =
        %|subset1| + |subset2|
    end
    %sample n times and count how many times a particular tree is picked
    sampled_tree = randsample(1:m, n, true);
    tree_counted = histc(sampled_tree, 1:m);
    
    %loop over all trees and sample points from the corresponding trees
    %with associate counts
    num_sampled = 0;
    
    fprintf('Resampling %d points...\n',n);
    %starting resampling
    for ctree = 1:m
        
        tree = sampler{ctree}; %The current tree for resampling
        prob = exp(sampler_prob{ctree}); %the probability for each leaf
        
        %Sample nodes with replacement tree_counted(ctree) times
        sampled_nodes = randsample(1:length(tree),tree_counted(ctree),true, prob);
        
        %Count how many times a unique node is picked
        unique_nodes = unique(sampled_nodes);
        nodes_count = histc(sampled_nodes, unique_nodes);
        
        %looping over all nodes picked to resample real samples
        for node_index = 1:length(unique_nodes)
            node = tree{unique_nodes(node_index)}; %the current sampled node
            if option.resample_data_only 
                %to resample using the observed data. data_sampled contains
                %rows of node.mark that have been sampled
                data_sampled = randsample(1:size(node.point,1),nodes_count(node_index), true);
                
                %Now we added the resampled points to y
                for k = data_sampled
                    data_loc = node.point(k,:);                
                    y(num_sampled + 1,:) = option.list(N(data_loc(2)) + data_loc(1),:);
                    num_sampled = num_sampled + 1;                    
                end
                
            elseif option.local_gaussian_smoothing && all(isfinite(node.cov(:)))
                %to resample using the normal density
                y(num_sampled+1:num_sampled+nodes_count(node_index),:) = ...
                    mvnrnd(node.mean,node.cov,nodes_count(node_index));
                num_sampled = num_sampled + nodes_count(node_index);
                
            else
                %to resample from the uniform distribution
                l = node.area(:,2) - node.area(:,1);
                y(num_sampled+1:num_sampled+nodes_count(node_index),:) =...
                    (ones(nodes_count(node_index),1)*l').*rand(...
                    nodes_count(node_index),d) + repmat(node.area(:,1)',nodes_count(node_index),1);
   
                num_sampled = num_sampled + nodes_count(node_index);
            end
        end
    end 
end