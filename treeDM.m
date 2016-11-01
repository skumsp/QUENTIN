function DM_tree = treeDM(tree)
n = size(tree,1);
DM_tree = zeros(n,n);
for i = 1:n
    child1 = tree(i,1);
    child2 = tree(i,2);
    if child1+child2 > 0
        DM_tree(i,child1) = tree(i,3);
        DM_tree(child1,i) = tree(i,3);
        DM_tree(i,child2) = tree(i,4);
        DM_tree(child2,i) = tree(i,4);
    end
end