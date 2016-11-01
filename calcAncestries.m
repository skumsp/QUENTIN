function Anc = getAncestries(TransNet)
nSamp = size(TransNet,2);
Anc = zeros(nSamp,nSamp);

G = digraph(TransNet);
for i = 1:nSamp
    v = bfsearch(G,i);
    for j=1:size(v,1)
        u = v(j);
        if u~=i
            Anc(i,u) = 1;
        end
    end
end