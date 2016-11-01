function TransNet = treeToTransNet(tree,DSamp)
nTreeVert = size(tree,1);
nVert = (nTreeVert + 1)/2;
TransNet = zeros(nVert,nVert);
for i=1:nTreeVert
    l = tree(i,3);
    child1 = tree(i,1);
    child2 = tree(i,2);
    if child1 + child2 == 0
        continue;
    end
    l1 = tree(child1,3);
    l2 = tree(child2,3);
    if l==l1
        TransNet(l,l2) = DSamp(l,l2);
    end
    if l == l2
        TransNet(l,l1) = DSamp(l,l1);
    end
end