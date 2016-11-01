function [DM_tree,vertLabels] = treeDMdirect(tree,names_pat)
n = size(tree,1);
DM_tree = zeros(n,n);
vertLabels = cell(1,n);
for i = 1:n
    l = tree(i,5);
    vertLabels{i} = ['n',int2str(i),'_',names_pat{l}];
    child1 = tree(i,1);
    child2 = tree(i,2);
    if child1+child2 > 0
        DM_tree(i,child1) = tree(i,3);
        DM_tree(i,child2) = tree(i,4);
    end
end