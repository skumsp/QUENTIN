function [newTree,w] = calcLabelsPhyloTree(tree,DM)
nVert = size(tree,1);
s = zeros(1,nVert-1);
t = zeros(1,nVert-1);
e = 1;
for i = 1:nVert
    if tree(i,1) > 0
        s(e) = i;
        t(e) = tree(i,1);
        e = e+1;
    end
    if tree(i,2) > 0
        s(e) = i;
        t(e) = tree(i,2);
        e = e+1;
    end
end
G = graph(s,t);
order = flipud(dfsearch(G,nVert));
outdeg = (tree(:,1)>0) + (tree(:,2)>0);
internal = outdeg > 0;
w = 1;
for i = 1:nVert
    v = order(i);
    if internal(v) == 1
        child1 = tree(v,1);
        child2 = tree(v,2);
        l1 = tree(child1,3);
        l2 = tree(child2,3);
        if DM(l1,l2) < DM(l2,l1)
            tree(v,3) = l1;
        else
            tree(v,3) = l2;
        end
        if (DM(l1,l2) == intmax) && (DM(l2,l1) == intmax)
            w = 0;
        end
    end
end
newTree = tree;
