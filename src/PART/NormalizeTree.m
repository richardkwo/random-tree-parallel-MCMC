function [tree, sampler, sampler_prob] = NormalizeTree(tree, MCdraws, log_total_sum)
    %NormalizeTree is a function that takes an unnormalized tree
    %and return a normalized tree. When NormalizeTree is called with
    %only two arguments, it will first search the whole tree, compute the
    %treeNode.log_prob (leaf probability) and treeNode.log_density (leaf density)
    %for each leaf with information in treeNode.point and update treeNode.mean 
    %and treeNode.var using MCdraws. It will simultaneously add up
    %the total probability and use the obtained total probability to
    %normalize the tree, and return the normalized tree + normalized
    %flattened tree(sampler) and the corresponding (sampler_prob)
    % 
    %When NormalizeTree is called with three arguments, it searches over
    %tree and normalize the tree but only return the normalized tree
    %without all computation as well as sampler and sampler_prob
    %
    %It is called by NormalizeForest to handle forest input
    %
    %Call:
    %tree = Normalize(tree, [], log_total_sum)
    %[tree, sampler, sampler_prob] = Normalize(tree, MCdraws)
    %
    %Arguments:
    %
    %tree: an output by buildTree, or an element of trees cells outputted by
    %   buildForest, OneStageMCMC or MultiStageMCMC, or any copied structure.
    %MCdraws: A list of subsets containing posterior samples on different
    %   subsets
    %log_total_sum: the logarithm of sums of exp(leaf.log_prob)
    %
    %Outputs:
    %
    %trees: Normalized random partition tree for density
    %   evaluation.
    %sampler: Normalized flattened tree for resampling (only outputted when
    %   the number of inputs is two)
    %sampler_prob: Corresponding probability for sampler (only outputted when
    %   the number of inputs is two)
    %
    %See also:
    %buildTree, NormalizeForest, buildForest
    %
    
    if nargin == 2 
        %Two inputs, need to compute all components for leaves and then
        %normalize the tree and output the sampler
        
        total_set = length(MCdraws);%total subsets
        
        %Initialize sampler, sampler_prob, dimension and subset posteiror sizes
        sampler = {}; 
        sampler_prob = [];
        d = size(MCdraws{1},2);
        % size of each subset
        N = zeros(1,total_set);
        for i = 1:total_set
           N(i) = size(MCdraws{i},1);
        end
        
        %Initialize the total probability
        %total_prob = 0;
        
        %We use an stack to write the non-recursion way of tree traverse
        stack = tree;
        
        while ~isempty(stack) 
            %when stack is not empty, pop up the last element
            current_node = stack(end);
            stack = stack(1:end-1); %pop the last one
            
            if isempty(current_node.left) && isempty(current_node.right)
                %If this element is leaf already, do the calculation...
                
                l = current_node.area(:,2) - current_node.area(:,1);%Size of the area
                
                %If multiple subsets
                if total_set > 1
                    %count numbers of points in different subsets.
                    if ~isempty(current_node.point)
                        c = histc(current_node.point(:,2), 1:total_set); % treeNode.point(i,j) = i-th sample in subset j
                    else
                        c = zeros(1, total_set);
                    end
                    counts = c + 0.01;%bound away from 0
                    
                    %compute the probability and the density
                    current_node.log_prob = sum(log(counts)) - sum(log(N)) - (total_set - 1)*sum(log(l)); %leave normalizing to next step in logarithm
                    current_node.log_density = current_node.log_prob - sum(log(l)); %divide by the area in logarithm
                
                %If only one subset
                else
                    n = size(current_node.point,1);
                    current_node.log_prob = log(n) - log(N);
                    current_node.log_density = current_node.log_prob - sum(log(l));
                end
                
                %Adding up the unnormalized probability
                %total_prob = total_prob + exp(current_node.prob);
                
                %Next, to compute the node.mean and node.var
                if total_set>1
                    %When multiple subsets. Because we just counted the
                    %number of points within each subset in "c". We can
                    %make use of that to pick the original data point from
                    %MCdraws.
                    
                    num_point = size(current_node.point,1); %total number of points
                    list = zeros(num_point,d); %temporal list to cached these points
                    cached = 0; %number of poins cached already
                    
                    %looping over all subsets
                    for cached_set = 1:total_set
                        if c(cached_set)>0
                            %If there are points from cached_set, we obtain
                            %index from node.point, transform it back to
                            %their id in the original subset
                            index = (cached+1): (cached + c(cached_set)); 
                            list(index,:) = MCdraws{cached_set}(current_node.point(index,1),:);
                            cached = cached + c(cached_set);
                        end
                    end
                    
                    %Compute mean and variance
                    current_node.mean = mean(list);
                    current_node.var = var(list)/total_set;
                else
                    %when there is only one subset
                    current_node.mean = mean(MCdraws{1}(current_node.point(:,1),:));
                    current_node.var = var(MCdraws{1}(current_node.point(:,1),:));
                end
                
                %argmenting the sampler and sampler probability
                sampler = [sampler, {current_node}];
                sampler_prob = [sampler_prob, current_node.log_prob];
                
            %If the current node is not a leaf, we push its left and right
            %child into the stack.
            else
                stack = [stack, current_node.left, current_node.right];
            end
        end
        
        %We normalize the sampler_prob by using total_probability
        max_prob = max(sampler_prob);
        
        %Now we call NormalizeTree again but using three inputs to get tree
        %normalized by log(total_prob)
        tree = NormalizeTree(tree, [], max_prob + log(sum(exp(sampler_prob - max_prob))));
        sampler_prob = sampler_prob - max_prob - log(sum(exp(sampler_prob - max_prob)));
                
    else
        %When the function is called with three inputs, we just use
        %log_total_sum to normalize it.
        
        %We use an iterative way to transverse the tree
        stack = tree;
        
        %Just to initialize to avoid any error information
        sampler = {};
        sampler_prob = [];
        
        %Transpassing the tree and update the node.prob and node.value
        while ~isempty(stack)
            if isempty(stack(end).left) && isempty(stack(end).right)
                % this is a leaf node
                stack(end).log_prob = stack(end).log_prob - log_total_sum;
                % 1-1 mapped: normalized density = normalized prob / (area of that block)
                stack(end).log_density = stack(end).log_density - log_total_sum;
                stack = stack(1:end-1); %pop the last one
            else
                stack = [stack(1:end-1), stack(end).left, stack(end).right];
            end
        end
    end
end