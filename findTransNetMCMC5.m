function [TransNets,TransNetTreesWeight,record,misc,stats] = findTransNetMCMC5(DSamp,nIter,maxdist,nEdgeModif,powers,intraHostCoeff,ndecr,initTree, debugMode,TransNetTrue,aux)

sigma = 0.1;
epsilon = 0.5;
delta = 0.0001;
gamma = 3;
npow = size(powers,2);
nEdgeModifCurr = nEdgeModif;
iEdgeModif = 1;
ipow = 1;
pow = powers(ipow);
stats = cell(1,2);

aux1 = [];
aux2 = [];
if ~isempty(aux)
    aux1 = aux(1,:);
    aux2 = aux(2,:);
end

nSamp = size(DSamp,1);
nVertTree = 2*nSamp-1;

% DSamp_dir = findDSampDir(DSamp,maxdist);
DSamp_dir = findDSampDir1(DSamp);

DSamp_symm= bsxfun(@min,DSamp,DSamp');


%tree - n x 3 matrix. First and second columns - children, third column -
%label

if isempty(initTree)

%     phyloTree = seqneighjoin(DSamp_symm);
    phyloTree = seqlinkage(DSamp_symm);
    % view(phyloTree)
    treeData = get(phyloTree);
    branchLists= treeData.Pointers;
    nodeNames = treeData.NodeNames;
    leafID = zeros(1,nSamp);
    for i=1:nSamp
        name = nodeNames{i};
        id = str2num(name(6:end));
        leafID(i) = id;
    end

    tree = zeros(nVertTree,3);
    for i=1:nSamp
        tree(i,3) = i;
    end
    for i = 1:(nSamp-1)
        if branchLists(i,1)<=nSamp
            tree(nSamp+i,1) = leafID(branchLists(i,1));
        else
            tree(nSamp+i,1) = branchLists(i,1);
        end
        if branchLists(i,2)<=nSamp
            tree(nSamp+i,2) = leafID(branchLists(i,2));
        else
            tree(nSamp+i,2) = branchLists(i,2);
        end
    end

    [tree,w] = calcLabelsPhyloTree(tree,DSamp);
end

% recordTree = tree;
TransNet = treeToTransNet(tree,DSamp);
k = estimateHubs(TransNet);
% k = 1;
% MatrToFit = TransNet;
MatrToFit = DSamp_dir;
AM_tree = treeAM(tree);
[obj1,treeWeight] = objTransNetPhyloFit3(AM_tree,tree,MatrToFit,intraHostCoeff,aux1);
% [obj1,treeWeight] = objTransNetPhyloFit4(AM_tree,tree,MatrToFit,intraHostCoeff,aux1);
obj2LinkNoPow = objTransNetQuadDeg2(TransNet,k,aux2);
obj2 = obj2LinkNoPow^pow;

% degDistr = getDegDistr(TransNet);
% obj2 = objTransNetDegDistrFit(degDistr,gamma);

record = w*obj1*obj2;
TransNetTrees = {tree};
TransNetTreesWeight = {treeWeight};

compToTrue = -1;
if ~isempty(TransNetTrue) 
    compToTrue = compareTransNets(TransNetTrue,TransNet);
end
    

solDynam = [obj1 obj2 record compToTrue];

treeLink = tree;
obj1Link = obj1;
w_link = w;
treeWeightLink = treeWeight;
compToTrueLink = compToTrue;

% names_pat = cellfun(@num2str, num2cell(1:nSamp), 'UniformOutput', false);
% bg = biograph(sparse(TransNet),names_pat,'ShowWeights','on');
% view(bg);

['start local search']

i = 0;
while i <= nIter*nEdgeModif
    i = i+1;
    [num2str(i) ' ' num2str(record) ' ' num2str(compToTrue)]
    
    if iEdgeModif == nIter+1
        iEdgeModif = 1;
%         nEdgeModif = nEdgeModif - 1;
        nEdgeModifCurr = nEdgeModifCurr - 1;
    end
    
    newtree = modifyTree1(tree,nEdgeModifCurr,nSamp,'NNI');
    if isempty(newtree)
        continue;
    end
    
    AM_newtree = treeAM(newtree);
    [newtree,w] = calcLabelsPhyloTree(newtree,DSamp);
    TransNetNew = treeToTransNet(newtree,DSamp);
%     MatrToFit = TransNetNew;
    MatrToFit = DSamp_dir;
    [obj1,newtreeWeight] = objTransNetPhyloFit3(AM_newtree,newtree,MatrToFit,intraHostCoeff,aux1);
%     [obj1,newtreeWeight] = objTransNetPhyloFit4(AM_newtree,newtree,MatrToFit,intraHostCoeff,aux1);
%     obj1 = objTransNetPhyloFit(AM_tree,D,sigma,epsilon);
    obj2 = objTransNetQuadDeg2(TransNetNew,k,aux2)^pow;
    
%     k_new = estimateHubs(TransNetNew);
%     degDistr = getDegDistr(TransNet);
%     obj2 = objTransNetDegDistrFit(degDistr,gamma);
    
    val = w*obj1*obj2;
    
%     if ~isempty(TransNetTrue)
%         compToTrue1 = compareTransNets(TransNetTrue,TransNetNew);
%     end
    
    prob = val/record;
%     p = rand;
    p = 1;
    if p <= prob
        tree = newtree;
    end
    if (val > record + delta) %&& (obj1 >= solDynam(end,1))
%         tree = newtree;
%         recordTree = newtree;
        record = val;
        TransNetTrees = {newtree};
        TransNetTreesWeight = {newtreeWeight};
        if ~isempty(TransNetTrue)
            compToTrue = compareTransNets(TransNetTrue,TransNetNew);
        end
        solDynam = [solDynam; obj1 obj2 obj1*obj2 compToTrue];
        
        if (ipow < npow) && (size(solDynam,1) >= ndecr)
            obj1prew = solDynam((end-ndecr+1):end,1);
            if all(diff(obj1prew)<0) 
   % decrease power and restart
                ipow = ipow + 1;
                pow = powers(ipow);
                nEdgeModifCurr = nEdgeModif;
                iEdgeModif = 1;
                i = 1;
                tree = treeLink;
                obj1 = obj1Link;
                obj2 = obj2LinkNoPow^pow;
                record = w_link*obj1*obj2;
                treeWeight = treeWeightLink;
                compToTrue = compToTrueLink;
                solDynam = [obj1 obj2 record compToTrue];
                TransNetTrees = {tree};
                TransNetTreesWeight = {treeWeight};
            end
        end
        
        
    end
    if abs(val-record) < delta
        iscont = false;
        nopt = size(TransNetTrees,2);
        for c = 1:nopt
            if isequal(TransNetTrees{c},newtree)
                iscont = true;
                break;
            end
        end
        if ~iscont
            TransNetTrees{nopt+1} = newtree;
            TransNetTreesWeight{nopt+1} = newtreeWeight;
        end
    end    
    iEdgeModif = iEdgeModif + 1;
end

solDynam
stats{1} = solDynam;
stats{2} = pow;
stats{3} = k;
misc = [mean(solDynam(:,1)) std(solDynam(:,1)); mean(solDynam(:,2)) std(solDynam(:,2))]

nopt = size(TransNetTrees,2);
TransNets = {};
for i = 1:nopt
    tree = TransNetTrees{i};
    AM_tree = treeToTransNet(tree,DSamp);
    isExist = 0;
    for j = 1:size(TransNets,2)
        if isequal(AM_tree,TransNets{j})
            isExist = 1;
            break;
        end
    end
    if isExist == 0
        TransNets{end+1} = AM_tree;
    end
end

% TransNet = treeToTransNet(recordTree,DSamp);