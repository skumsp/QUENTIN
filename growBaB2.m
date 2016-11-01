function [] = growBaB2(V_arb_label,E_arb,F,DM)
% global variables: Arbs
global Arbs;
global bestArbDM;
global record;
global gamma;
global arbCount;
global bound;
global s;
global disconCount;

nVert = size(DM,2);

if max(V_arb_label) > bound
%     ['branch cutted by dist']
    return;
end

% if size(E_arb,2) > 0
%     w = 0;
%     for i=1:size(E_arb,2)
%         e = E_arb{i};
%         w = w+e(3);
%     end
%     if w > bound
% %     ['branch cutted by dist']
%         return;
%     end
% end


% lb = getLowerBound(E_arb,DM,gamma);
% if lb >= record
% %         ['branch cutted by distr']
%         return;
% end

% G = digraph(DM);
% v = bfsearch(G,s);
% if size(v,1) < nVert
% %     ['branch cutted by conn']
%     disconCount = disconCount + 1;
%     return;
% end

V_arb = V_arb_label >= 0;
if sum(V_arb,2) == nVert
    arbCount = arbCount + 1;
    ['found ', int2str(arbCount)]
%     if arbCount == 5000
%         elapsedTime = toc;
%         [int2str(elapsedTime)]
%     end
    
%     Arbs{end+1} = E_arb; 
    DM = eListToAM(E_arb,nVert);
    degDistr = getDegDistr(DM);
    val = objTransNet(degDistr,gamma);
    if val < record
        record = val;
        bestArbDM = DM;
    end
    return;
end
% for i = 1:size(F,1)
    if size(F,1) == 0
        return;
    end
    e = F(1,:);
    u = e(1);
    v = e(2);
    E_arb{end+1} = e;
    V_arb_label(v) = V_arb_label(u) + e(3);
    
    mFnew = size(F,1);
    removeEdges = zeros(size(F,1),1);
    for i=1:size(F,1)
        e1 = F(i,:);
        w = e1(1);
        if (e1(2) == v) && (V_arb_label(w) >= 0)
            removeEdges(i) = 1;
            mFnew = mFnew-1;
        end
    end
    for w=1:nVert
        if (DM(v,w) > 0) && (V_arb_label(w) < 0)
            mFnew = mFnew + 1;
        end
    end
   
    F_new = zeros(mFnew,3);
    j = 1;
    for i=1:size(F,1)
        if removeEdges(i) == 0
            F_new(j,:) =  F(i,:);
            j = j+1;
        end
    end
    for w=1:nVert
        if (DM(v,w) > 0) && (V_arb_label(w) < 0)
            F_new(j,:) = [v w DM(v,w)];
            j = j+1;
        end
    end

    
    growBaB2(V_arb_label,E_arb,F_new,DM);
    E_arb(end) = [];
    V_arb_label(v) = -1;
    DM(u,v) = 0;
    
%     G = digraph(DM);
%     v = bfsearch(G,s);
%     if size(v,1) == nVert
    isConn = isConnected(V_arb_label,DM);
    if isConn
        F_new = F(2:end,:);
        growBaB2(V_arb_label,E_arb,F_new,DM); 
%     else
%         disconCount = disconCount + 1;
    end
%     ['Leaf']
% end


