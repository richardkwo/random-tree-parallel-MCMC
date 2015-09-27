function y = treeDensity(x, trees, dens)
    %treeDensity is a function taking output from buildForest,
    %OneStageMCMC and MultiStageMCMC to evaluate the combined density for a given
    %point.
    %
    %Call:
    %y = treeDensity(x, trees)
    %y = treeDensity(x, trees, dens)
    %
    %Arguments:
    %
    %x: is a N x d matrix containing all points that density value need be
    %   evaluated. N is number of points and d is the dimension.
    %trees: forest of random partition tree (or flattened tree, depending 
    %   on "dens") outputed ffrom buildForest, OneStageMCMC and MultiStageMCMC.
    %   When dens is "uniform", trees are
    %   "trees" outputted from the three functions. When dens is "normal",
    %   trees are "sampler" outputted from the three functions.
    %dens: choices of different density functions on leaves. Options are
    %"uniform" and "normal". Default is "uniform".
    %
    %Outputs:
    %
    %y: the density values at x
    %
    %See also:
    %buildForest, OneStageMCMC, MultiStageMCMC, treeSampling
    %
    
    
    m = length(trees); %number of tree
    [n, d] = size(x); %number of points and the dimension
    
    %Setting the default value for parameters
    if nargin <3
        dens = 'uniform';
    end
    
    %If the choice is uniform density on each leaf
    if strcmp(dens, 'uniform')
        %setting up waitbar...
        h = waitbar(0,'Density evaluation in parallel (uniform) ...');
        
        %Initializing z
        z = zeros(m,n);
        
        % parallel loop over trees
        for i=1:m
            for k=1:n
                Node = trees{i};
                while ~isnan(Node.dim)
                    if x(k, Node.dim) <= Node.value
                        Node = Node.left;
                    else
                        Node = Node.right;
                    end
                end
                %Adding the density to y(k)
                z(i, k) = exp(Node.value);
            end
        end
        
        %Averaging over all trees
        y = mean(z, 1);

        %updating waitbar...
        waitbar(1, h, 'Evaluating densities ...');
        close(h);
    end
    
    %If the choice is to use normal density on each leaf
    if strcmp(dens,'normal')
        %setting up the waitbar
        h = waitbar(0,'Starting density evaluation (normal fit) ...');
                
        %looping over all points
        for t = 1:n
        %looping over all trees
        for i = 1:m
            y1 = 0;
            %looping over all leaves in that tree
            for k = 1:length(trees{i})
                %l = trees{i}{k}.area(:,2) - trees{i}{k}.area(:,1);
                y1 = y1 + mvnpdf(x(t,:),trees{i}{k}.mean',trees{i}{k}.cov);
            end
            y1 = y1/length(trees{i});
            y(t) = y(t) + y1;
        end
        y(t) = y(t)/m;
        
        %updating waitbar
        waitbar(t/n, h, 'Evaluating densities (guassian fit) ...');
        end
        close(h);
    end
end