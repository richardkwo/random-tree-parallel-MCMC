function plot_tree_blocks(tree, sub_chain, varargin)
% visualize the blocks estimated by a tree by drawing pairwise block
% structures
% Call
% plot_tree_blocks(treeNode)
% plot_tree_blocks(tree_sampler, [])
% plot_tree_blocks(tree_sampler, sub_chain)
% plot_tree_blocks(tree_sampler, sub_chain, [p1 ... pn]), if only a few dimensions are plotted

if ~isempty(sub_chain)
    M = length(sub_chain);
    cc = hsv(M);
else
    cc = flipud(gray(1000));
end

if isa(tree, 'treeNode')
    % get its leaf nodes
else
    max_logdens = 0;
    n_leaf = length(tree);
    for t=1:n_leaf 
        if tree{t}.log_density>max_logdens
            max_logdens = tree{t}.log_density;
        end
    end

    assert(isa(tree{1}, 'treeNode'), 'wrong input');
    fprintf('Plotting %d blocks...\n', n_leaf);
    % tree is a set of leaf nodes (samplers)
    p = size(tree{1}.area, 1);
    if ~isempty(varargin)
        dims = varargin{1};
        assert(all(dims>=1) && all(dims<=p));
    else
        dims = 1:p;
    end
    a = length(dims);
    % pairwise plots
    for i=dims
        for j=dims
            hs = subplot(a,a,(i-1)*a+j);
            hold on;
            if i==j
                text(0.5, 0.5, ['parameter ',num2str(j)], 'Parent', hs);
            elseif i<j
                % iterate over and draw blocks
                for t=1:n_leaf
                    x = tree{t}.area(i,1);
                    y = tree{t}.area(j,1);
                    w = tree{t}.area(i,2) - tree{t}.area(i,1);
                    h = tree{t}.area(j,2) - tree{t}.area(j,1);
                    
                    % if sub_chains supplied, draw the points and rectangle
                    % boundary
                    if ~isempty(sub_chain)
                        rectangle('Position',[x y w h], 'EdgeColor','red');
                        points = zeros(size(tree{t}.point, 1),2);
                        for k=1:size(points,1)
                            c = tree{t}.point(k,2);
                            internal_index = tree{t}.point(k,1);
                            points(k,1) = sub_chain{c}(internal_index, i);
                            points(k,2) = sub_chain{c}(internal_index, j);
                        end
                        plot(points(:,1), points(:,2), '.', 'color', cc(c,:), 'LineWidth', 0.1);
                    else
                        % draw filled rectangle
%                         dens = (1+tanh((tree{t}.log_prob - max_logdens)/50))/2
%                         rectangle('Position',[x y w h],...
%                             'FaceColor', cc(floor(dens * size(cc,1))+1,:));
                        rectangle('Position',[x y w h], 'LineWidth', 0.1);
                    end
                end
            end
            hold off;
            xlabel(['parameter ', num2str(i)]);
            ylabel(['parameter ', num2str(j)]);
        end
    end
end
    
end