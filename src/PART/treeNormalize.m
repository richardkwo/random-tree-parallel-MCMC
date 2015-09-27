function [sum2, prob] = treeNormalize(tree, sum1)
    if nargin == 1
        %to get the sum and the probability vector
        if isnan(tree.left) && isnan(tree.right)
            sum2 = tree.value;
            prob = [tree.prob];
        else
            [sum2_left, prob_left] = treeNormalize(tree.left);
            [sum2_right, prob_right] = treeNormalize(tree.right);
            sum2 = sum2_left + sum2_right;
            prob = [prob_left prob_right];
        end
        
    elseif nargin == 2
        %to normalize using sum1
        sum2 = 0;
        prob = 0;
        if isnan(tree.left) && isnan(tree.right)
            tree.value = tree.value/sum1;
        else
            treeNormalize(tree.left, sum1);
            treeNormalize(tree.right, sum1);
        end
    end
end