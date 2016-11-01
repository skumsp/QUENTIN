function [val,tree_weight] = objTransNetPhyloFit4(AM_tree,tree,DM,intraHostCoeff,aux)
%     sigma = 0.1;
%     epsilon = 0.5;

%     interHostCoeff = 0.85;
        n = size(DM,2);
        [v1, v2] = find(triu(AM_tree));
        E = [v1,v2];
        nEdges = size(E,1);
        
        intra_edges = zeros(1,nEdges);
        nIntraEdges = 0;
        for e = 1:nEdges
                edge = E(e,:);
                if tree(edge(1),3) == tree(edge(2),3)
                    intra_edges(e) = 1;
                    nIntraEdges = nIntraEdges + 1;
                end
        end
        E_intra = E(find(intra_edges > 0));
        
        nPairs = size(find(DM > 0),1);
        C = zeros(nPairs,nEdges);
        d = zeros(nPairs,1);

        pathsPairs = {};
        p = 0;
        for i=1:n
            for j = 1:n
                if DM(i,j) > 0
                    p = p + 1;
    %                 d(p) = DM(i,j);
                    d(p) = 1;
                    [dist, path, pred] = graphshortestpath(sparse(AM_tree), i, j);
                    pathsPairs{p} = path;
                    for v = 1:(size(path,2)-1)
                        if path(v) < path(v+1)
                            edge = [path(v) path(v+1)];
                        else
                            edge = [path(v+1) path(v)];
                        end
                        e = find(ismember(E,edge,'rows'));
    %                     C(p,e) = 1;
                        C(p,e) = 1/DM(i,j);
                    end
                end
            end
        end
        C = [C zeros(nPairs,nIntraEdges)];

        deg = sum(AM_tree,1);
        root = find(deg == 2);
        leafs = find(deg == 1);
        nleafs = size(leafs,2);
        [dists, pathsRoot, preds] = graphshortestpath(sparse(AM_tree), root);
        pathsRoot = pathsRoot(leafs);
        charVectPaths = zeros(nleafs,nEdges);
        for i = 1:nleafs
             path = pathsRoot{i};
             for v = 1:(size(path,2)-1)
                        if path(v) < path(v+1)
                            edge = [path(v) path(v+1)];
                        else
                            edge = [path(v+1) path(v)];
                        end
                        e = find(ismember(E,edge,'rows'));
                        charVectPaths(i,e) = 1;
             end
        end
        charVectPathsIntra = charVectPaths.*repmat(intra_edges,nleafs,1);
        charVectPaths = [charVectPaths charVectPathsIntra(:,E_intra)];
        
        Aineq = zeros(nIntraEdges,nEdges  + nIntraEdges);
        bineq = zeros(nIntraEdges,1);
        for i = 1:nIntraEdges
            e = E_intra(i);
            Aineq(i,e) = -intraHostCoeff;
            Aineq(i,nEdges + i) = 1;
        end

        Aeq = zeros(nleafs-1,nEdges + nIntraEdges);
        beq = zeros(nleafs-1,1);
        for i=1:(nleafs-1)
            Aeq(i,:) = charVectPaths(i,:)-charVectPaths(i+1,:);
        end
        lb = zeros(nEdges+nIntraEdges,1);
        ub = Inf*ones(nEdges+nIntraEdges,1);

%         x = lsqlin(C,d,[],[],Aeq,beq,lb,ub);
        x = lsqlin(C,d,Aineq,bineq,Aeq,beq,lb,ub);


        tree_weight = zeros(size(tree,1),5);
        for e=1:nEdges
            u = E(e,1);
            v = E(e,2);
            if (tree(u,1) == v) | (tree(u,2) == v)
                par = u;
                child = v;
            else
                par = v;
                child = u;
            end
            if tree_weight(par,1) == 0
                tree_weight(par,1) = child;
                tree_weight(par,3) = x(e);
            else
                tree_weight(par,2) = child;
                tree_weight(par,4) = x(e);
            end        
        end
        for i = 1:size(tree,1)
            tree_weight(i,5) = tree(i,3);
        end

        DSamp_tree = zeros(n,n);
        p = 0;
        for i=1:n
            for j = 1:n
                if DM(i,j) > 0
                    p = p + 1;
                    pathEdges = find(C(p,1:nEdges) > 0);
                    DSamp_tree(i,j) = sum(x(pathEdges),1);
                end
            end
        end

        vec1 = DM(:);
        vec2 = DSamp_tree(:);
        ind = find(vec1 > 0);
        vec1 = vec1(ind);
        vec2 = vec2(ind);
        p = size(ind,1);
        val = max(corr(vec1,vec2),0);
        
%     val = corr(vec1,vec2,'type','Kendall');
         
    

%     val = 1;
%     p = 0;
%     for i=1:n
%         for j = 1:n
%             if DM(i,j) > 0
%                 p = p + 1;
%                 path = pathsPairs{p};
%                 mu_p = 0;
%                 for v = 1:(size(path,2)-1)
%                     if path(v) < path(v+1)
%                         edge = [path(v) path(v+1)];
%                     else
%                         edge = [path(v+1) path(v)];
%                     end
%                     e = find(ismember(E,edge,'rows'));
%                     mu_p = mu_p + x(e);
%                 end
%                 sigma_p = sigma*(size(path,2)-1);
%                 prob = normcdf([DM(i,j)-epsilon, DM(i,j)+epsilon],mu_p,sigma_p);
%                 val = val*(prob(2)-prob(1));
%             end
%         end
%     end

%     val
% 
%     DM_tree = zeros(2*n-1,2*n-1);
%     DM_tree_time = zeros(2*n-1,2*n-1);
%     intra_edges_set = find(intra_edges > 0);
%     for i=1:(2*n-1)
%         for j=(i+1):(2*n-1)
%             if AM_tree(i,j) == 1
%                 edge = [i j];
%                 e = find(ismember(E,edge,'rows'));
%                 DM_tree(i,j) = x(e);
%                 DM_tree_time(i,j) = x(e);
%                 if intra_edges(e) == 1
%                     e1 = find(intra_edges_set == e);
%                     DM_tree_time(i,j) = DM_tree_time(i,j)+x(nEdges + e1);
%                 end
%             end
%         end
%     end
%     bgTree = biograph(DM_tree',[],'ShowWeights','on');
%     view(bgTree);
%     
%     bgTree = biograph(DM_tree_time',[],'ShowWeights','on');
%     view(bgTree);
%     
%     ['finished']