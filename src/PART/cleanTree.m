function [trees, numNode] = cleanTree(trees)
    %This is a function that cleans a given forest of random partition
    %tree. In particular, it cleans up node.log_prob, node.log_density and node.point
    %but leaves all gates(partition point) and area information unchanged
    %(i.e., the tree structure is preserved)
    
    m = length(trees);
    h = waitbar(0,'Starting cleaning trees...');
    for i = 1:m
        Node = trees{i};
        stack = Node;
        while ~isempty(stack)
            current_node = stack(end);
            stack = stack(1:end-1); %pop the last one
            if isempty(current_node.left) && isempty(current_node.right)
                current_node.log_prob = nan;
                current_node.log_density = nan;
                current_node.point = [];
                current_node.mean = nan;
                current_node.cov = nan;
            else
                stack = [stack, current_node.left, current_node.right];
            end
        end
        waitbar(i/m, h, ['Cleaning', num2str(i), 'th tree...']);
    end
    close(h);
end