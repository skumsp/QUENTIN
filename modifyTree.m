function newtree = modifyTree(tree,steps,nSamp,method)
    u = randi(nSamp-1)+nSamp;
    j = randi(2);
    v = tree(u,j);
    newtree = tree;
    if j == 1
        k = randi(2);
        w = newtree(u,2);
        s = newtree(v,k);
        if newtree(v,k) == 0
            newtree = [];
            return;
        end
        newtree(u,2) = newtree(v,k);
        newtree(v,k) = w;
    end
    if j == 2
        k = randi(2);
        w = newtree(u,1);
        if newtree(v,k) == 0
            newtree = [];
            return;
        end
        newtree(u,1) = newtree(v,k);
        newtree(v,k) = w;
    end