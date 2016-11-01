function [DSamp, TransNet,TransTree,iter, AMCont] = simOutbreak(nSamp,alpha, infRate,maxIter,collectAfterTrans,graphAlg, debugMode)

DSamp = zeros(nSamp,nSamp);
TransNet = zeros(2*nSamp-1,3);

% infRate = 5/100; % the probability of infection transfer per link per unit time. 

if strcmp(graphAlg,'muchnik')
    Alpha = -alpha;
    % XAxis = [1 : nSamp-1];
    % YAxis = [1 : nSamp-1].^(Alpha+1);
    XAxis  = unique(round(logspace(0,log10(nSamp),25)));
    YAxis  = unique(round(logspace(0,log10(nSamp),25))).^(Alpha+1);
    Graph = mexGraphCreateRandomGraph(nSamp,XAxis,YAxis,0); 

    % Degrees = GraphCountNodesDegree(Graph);
    % h1 = figure;
    % % incoming:
    % [y x] = hist(Degrees(:,2),unique(Degrees(:,2)));
    % loglog(x,y/sum(y),'*r');
    % hold on
    % % outgoing
    % [y x] = hist(Degrees(:,3),unique(Degrees(:,3)));
    % loglog(x,y/sum(y),'dg');
    % % expected distribution:
    % loglog(XAxis,YAxis/sum(YAxis),':b');
    % xlabel('k,Degree');
    % ylabel('P(k)');
    % title('Node Degree Distribution');
    % legend({'Incoming','Outgoing','Expected'});

    nSamp = GraphCountNumberOfNodes(Graph);
    Econt_dir = Graph.Data;
    Econt = zeros(size(Econt_dir,1)/2,2);
    e1 = 1;
    for e = 1:size(Econt_dir,1)
        if Econt_dir(e,1) < Econt_dir(e,2)
            Econt(e1,:) = Econt_dir(e,1:2);
            e1 = e1 + 1;
        end
    end

    AMCont = eListToAM_matr(Econt,nSamp,true);
end

if strcmp(graphAlg,'octave')
    Econt = preferentialAttachment (nSamp, 1);
    AMCont = edgeL2adj(Econt);

%     AMCont = symmetrize(PriceModel(nSamp));

    nSamp = size(AMCont,1);
    
%     bg = biograph(sparse(triu(AMCont)),[],'ShowArrows','off');
%     view(bg);
%     
    Econt = buildSmaxGraph(degrees(AMCont));
    AMCont = symmetrize(edgeL2adj(Econt));
%     
%     
%     bg = biograph(sparse(triu(AMContmax)),[],'ShowArrows','off');
%     view(bg);
end

perm = randperm(nSamp);
AMCont = AMCont(perm,perm);
[S, C] = graphconncomp(sparse(AMCont), 'Directed', false);
maxcompsize = 0;
maxcomp = 0;
for i = 1:S
    compsize = sum(C == i,2);
    if compsize > maxcompsize
        maxcompsize = compsize;
        maxcomp = i;
    end
end
comp = find(C == maxcomp);
AMCont = AMCont(comp,comp);
nSamp = maxcompsize;
nEdges = sum(sum(AMCont,1),2)/2;

Econt = zeros(nEdges,2);
e1 = 1;
for i=1:nSamp
    for j = (i+1):nSamp
        if AMCont(i,j) == 1
            Econt(e1,1) = i;
            Econt(e1,2) = j;
            e1 = e1 + 1;
        end
    end
end
nEdges = size(Econt,1);


infected = zeros(1,nSamp);
infTime = zeros(1,nSamp);
e = Econt(randi(nEdges),:);
infected(e(1)) = 1;
infected(e(2)) = 1;
infTime(e(1)) = 0;
infTime(e(2)) = 1;

%child child st st fin fin label
treeSF = [2 3 1 1 0 0 e(1); 0 0 0 0 0 0 e(1); 0 0 0 0 0 0 e(2)];
parents = [0 1 1];

iter = 2;
nInf = 2;

while ((iter <= maxIter)&&(nInf < nSamp))
%         infSamp = find(infected == 1);
%         u = infSamp(randi(size(infSamp,2)));
%         neigh = adjLists{u};
%         v = neigh(randi(size(neigh,2)));
        [iter nInf]
        e = Econt(randi(nEdges),:);
        u = 0;
        v = 0;
        if (infected(e(1)) == 1) && (infected(e(2)) == 0)
            u = e(1);
            v = e(2);
        end
        if (infected(e(2)) == 1) && (infected(e(1)) == 0)
            u = e(2);
            v = e(1);
        end 
        if u == 0
            iter = iter + 1;
            continue;
        end
        p = rand;
        if p <= infRate
            infected(v) = 1;
            infTime(v) = iter;
            nInf = nInf + 1;
            treeu = find((treeSF(:,1) == 0)&(treeSF(:,7) == u));
            nTreeVert = size(treeSF,1);
            
            treeSF(treeu,1) = nTreeVert+1;
            treeSF(treeu,2) = nTreeVert+2; 
            treeSF(treeu,3) = iter;
            treeSF(treeu,4) = iter;
            treeSF(treeu,5) = 0;
            treeSF(treeu,6) = 0;
            par = parents(treeu);
            if treeSF(par,1) == treeu
                treeSF(par,5) = iter;
            end
            if treeSF(par,2) == treeu
                treeSF(par,6) = iter;
            end
            
            treeSF(nTreeVert+1,1) = 0;
            treeSF(nTreeVert+1,2) = 0; 
            treeSF(nTreeVert+1,3) = 0;
            treeSF(nTreeVert+1,4) = 0;
            treeSF(nTreeVert+1,5) = 0;
            treeSF(nTreeVert+1,6) = 0;
            treeSF(nTreeVert+1,7) = u;
            
            treeSF(nTreeVert+2,1) = 0;
            treeSF(nTreeVert+2,2) = 0; 
            treeSF(nTreeVert+2,3) = 0;
            treeSF(nTreeVert+2,4) = 0;
            treeSF(nTreeVert+2,5) = 0;
            treeSF(nTreeVert+2,6) = 0;
            treeSF(nTreeVert+2,7) = v;  
            
            parents(nTreeVert+1) = treeu;
            parents(nTreeVert+2) = treeu;
        end
        iter = iter + 1;
end
[iter nInf]
leafs = find(treeSF(:,1) + treeSF(:,2) == 0);
for i=1:size(leafs,1)
    l = leafs(i);
    par = parents(l);
    j = find(treeSF(par,1:2) == l);
    if j == 1
        treeSF(par,5) = iter + collectAfterTrans;
    end
    if j == 2
        treeSF(par,6) = iter + collectAfterTrans;
    end
end

% child child len len label
tree = zeros(size(treeSF,1),5);
tree(:,1:2) = treeSF(:,1:2);
tree(:,3) = treeSF(:,5) - treeSF(:,3);
tree(:,4) = treeSF(:,6) - treeSF(:,4);
tree(:,5) = treeSF(:,7);

infVert = find(infected > 0);
infTime = infTime(infVert);
for i = 1:size(tree,1)
    j = find(infVert == tree(i,5));
    tree(i,5) = j;
end
nSamp = nInf;

DSamp = -ones(nSamp,nSamp);
names = cellfun(@num2str, num2cell(1:nSamp), 'UniformOutput', false)
DM_tree = treeDM(tree);

[DM_tree_direct,labels] = treeDMdirect(tree,names);
if debugMode == 1
    bg = biograph(sparse(DM_tree_direct),labels,'ShowWeights','on');
    view(bg);
end

for i = 1:nSamp
    for j = (i+1):nSamp
        u = find((tree(:,1) == 0)&(tree(:,5)==i));
        v = find((tree(:,1) == 0)&(tree(:,5)==j));
        if (isempty(u)) || (isempty(v))
            continue;
        end
        [dist, path, pred] = graphshortestpath(sparse(DM_tree), u, v,'Directed', false);
        if infTime(i) < infTime(j)
            DSamp(i,j) = dist;
        else
            DSamp(j,i) = dist;
        end
    end
end

for i = 1:nSamp
    for j = 1:nSamp
        if DSamp(i,j) == -1
            DSamp(i,j) = intmax;
        end
        if i==j
            DSamp(i,j) = 0;
        end
    end
end

TransTree = tree(:,[1 2 5]);
TransNet = treeToTransNet(tree(:,[1 2 5]),DSamp);
['Done']



