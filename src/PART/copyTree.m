function newTrees = copyTree(trees)
    %This is a deepCopy function that copies a forest of random partition 
    %tree structure
    
    m = length(trees);
    newTrees = cell(1,m);
    h = waitbar(0,'Starting copying trees...');
    for i = 1:m
        Node = trees{i};
        newNode = copyNode(Node);
        newTrees{i} = newNode;
        stack = Node;
        newStack = newNode;
        while ~isempty(stack)
            current_node = stack(end);
            current_newNode = newStack(end);
            stack = stack(1:end-1); %pop the last one
            newStack = newStack(1:end-1);
            if ~isempty(current_node.left) || ~isempty(current_node.right)
                stack = [stack, current_node.left, current_node.right];
                newLeft = copyNode(current_node.left);
                newRight = copyNode(current_node.right);
                current_newNode.left = newLeft;
                current_newNode.right = newRight;
                newStack = [newStack, newLeft, newRight];
            end
        end
        waitbar(i/m, h, ['Copying', num2str(i), 'th tree...']);
    end
    close(h);
end

function new = copyNode(this)
    new = feval(class(this));
    p = properties(this);
    for i = 1:length(p)
        new.(p{i}) = this.(p{i});
    end
end