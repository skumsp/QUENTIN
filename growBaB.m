function [] = growBaB(V_arb,E_arb,F,AM)
% global variables: Arbs
global Arbs;
global bestArbAM;
global record;
global gamma;
global arbCount;
global timeBound;
nVert = size(AM,2);
% V_arb = V_arb_label >= 0;
if sum(V_arb,2) == nVert
    arbCount = arbCount + 1;
    ['found ', int2str(arbCount)]
    Arbs{end+1} = E_arb; 
    AM = eListToAM(E_arb,nVert) > 0;
    degDistr = getDegDistr(AM);
    val = objTransNet(degDistr,gamma);
    if val < record
        record = val;
        bestArbAM = AM;
    end
    return;
end
% for i = 1:size(F,1)
    if size(F,1) == 0
        return;
    end
    i = 1;
    e = F(i,:);
    F_new = F;
    F_new(i,:) = [];
    u = e(1);
    v = e(2);
    E_arb{end+1} = e;
    V_arb(v) = 1;
    for w=1:nVert
        if (AM(v,w) == 1) && (V_arb(w) == 0)
            F_new = [F_new; [v w]];
        end
        if (AM(w,v) == 1) && (V_arb(w) == 1) && (size(F_new,1) > 0)
            toRemoveEdges = [];
            for k=1:size(F_new,1)
                if (F_new(k,1) == w) && (F_new(k,2) == v)
                    toRemoveEdges = [toRemoveEdges; k];
                end
            end
            F_new(toRemoveEdges,:) = [];
        end
    end
    growBaB(V_arb,E_arb,F_new,AM);
    E_arb(end) = [];
    V_arb(v) = 0;
    AM(u,v) = 0;
    F_new = F;
    F_new(i,:) = [];
    growBaB(V_arb,E_arb,F_new,AM);  
%     ['Leaf']
% end


