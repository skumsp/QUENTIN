function [] = growBaB1(V_arb_label,E_arb,F,DM)
% global variables: Arbs
global Arbs;
global bestArbDM;
global record;
global gamma;
global arbCount;
global bound;
global s;

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

G = digraph(DM);
v = bfsearch(G,s);
if size(v,1) < nVert
%     ['branch cutted by conn']
    return;
end

% lb = getLowerBound(E_arb,DM,gamma);
% if lb >= record
% %         ['branch cutted by distr']
%         return;
% end


V_arb = V_arb_label >= 0;
if sum(V_arb,2) == nVert
    arbCount = arbCount + 1;
    ['found ', int2str(arbCount)]
    if arbCount == 5000
        elapsedTime = toc;
        [int2str(elapsedTime)]
    end
    
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
    F_new = F(2:end,:);
    u = e(1);
    v = e(2);
    E_arb{end+1} = e;
    V_arb_label(v) = V_arb_label(u) + e(3);
    for w=1:nVert
        if (DM(v,w) > 0) && (V_arb_label(w) < 0)
            F_new = [F_new; [v w DM(v,w)]];
        end
        if (DM(w,v) > 0) && (V_arb_label(w) >= 0) && (size(F_new,1) > 0)
            toRemoveEdges = [];
            for k=1:size(F_new,1)
                if (F_new(k,1) == w) && (F_new(k,2) == v)
                    toRemoveEdges = [toRemoveEdges; k];
                end
            end
            F_new(toRemoveEdges,:) = [];
        end
    end
    growBaB1(V_arb_label,E_arb,F_new,DM);
    E_arb(end) = [];
    V_arb_label(v) = -1;
    DM(u,v) = 0;
    F_new = F(2:end,:);
    growBaB1(V_arb_label,E_arb,F_new,DM);  
%     ['Leaf']
% end


