function newtree = modifyTree1(tree,steps,nSamp,method)
newtree = tree;
for i = 1:steps
    u = randi(nSamp-1)+nSamp;
    j = randi(2);
    v = newtree(u,j);
    newtree1 = newtree;
    if j == 1
        k = randi(2);
        w = newtree1(u,2);
        s = newtree1(v,k);
        if newtree1(v,k) == 0
            continue;
        end
        newtree1(u,2) = newtree1(v,k);
        newtree1(v,k) = w;
    end
    if j == 2
        k = randi(2);
        w = newtree1(u,1);
        if newtree1(v,k) == 0
            continue;
        end
        newtree1(u,1) = newtree1(v,k);
        newtree1(v,k) = w;
    end
    
    newtree = newtree1;
end