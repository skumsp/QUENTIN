function [] = grow(V_arb,E_arb,F,AM)
% global variables: Arbs
global Arbs;
global bestArbAM;
global record;
global gamma;
global arbCount;
nVert = size(AM,2);
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
b = false;
% import java.util.Stack;
% FF = java.util.Stack();
while ~b
%     if size(F,1) == 0
%         []
%     end
%     e = F.pop();
    e = F(end,:);
    F(end,:) = [];
    F_new = F;
    u = e(1);
    v = e(2);
    E_arb{end+1} = e;
    V_arb(v) = 1;
    for w=1:nVert
        if (AM(v,w) == 1) && (V_arb(w) == 0)
%             F_new.push([v w]);
            F_new = [F_new; [v w]];
        end
        if (AM(w,v) == 1) && (V_arb(w) == 1) && (size(F_new,1) > 0)
%             F_new.remove([w v]);
            toRemoveEdges = [];
            for k=1:size(F_new,1)
                if (F_new(k,1) == w) && (F_new(k,2) == v)
                    toRemoveEdges = [toRemoveEdges; k];
                end
            end
            F_new(toRemoveEdges,:) = [];
        end
    end
    grow(V_arb,E_arb,F_new,AM);
    E_arb(end) = [];
    V_arb(v) = 0;
    AM(u,v) = 0;
%     FF.push(e);
    AM_L = eListToAM(Arbs{end},nVert);
    [disc, pred, closed] = graphtraverse(sparse(AM_L),v,'Method','BFS');
    b = true;
    for w = 1:nVert
        if isnan(pred(w)) && (AM(w,v) == 1)
            b = false;
            break;
        end
    end
end
% while ~FF.isEmpty()
%     
% end


