function AM_tree = treeAM(tree)
n = size(tree,1);
AM_tree = zeros(n,n);
for i = 1:n
    child1 = tree(i,1);
    child2 = tree(i,2);
    if child1+child2 > 0
        AM_tree(i,child1) = 1;
        AM_tree(child1,i) = 1;
        AM_tree(i,child2) = 1;
        AM_tree(child2,i) = 1;
    end
end