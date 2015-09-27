classdef treeNode < handle
    % treeNode is either a gate (non-leaf nodes that splits into two), or a
    % node (leaf nodes that represents an equal-density block)
    properties
        dim = nan %The cutting dimension (for gates)
        value = nan %The cutting point (for gates)
        log_density = nan % The density (for nodes) (in log-scale)
        log_prob = nan %The node probability (for nodes) (in log-scale)
        point = [] %Indexes for original posterior samples. (i,j) means the ith point in jth subset
        area = nan %Area information
        cov = nan %variance of normal density approximation
        left = [] %reference to the left child
        right = [] %reference to the right child
        mean = nan %mean of normal density approximation
    end
    methods
        function new = copyNode(this)
            new = feval(class(this));
            p = properties(this);
            for i = 1:length(p)
                new.(p{i}) = this.(p{i});
            end
        end
    end
end