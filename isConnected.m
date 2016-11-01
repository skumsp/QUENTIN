function isConn = isConnected(V_arb,DM)
nVert = size(DM,1);
labels = V_arb >= 0;
visited = find(labels == 1);
nVisited = size(visited,2);
queue = zeros(1,nVert);
queue(1:nVisited) = visited;
head = 1;
tail = nVisited;
while head <= tail
    currvert = queue(head);
    for i = 1:nVert
        if (DM(currvert,i) > 0) && (labels(i) == 0)
            labels(i) = 1;
            tail = tail + 1;
            queue(tail) = i;
        end
    end
    head = head + 1;
end
if sum(labels,2) == nVert
    isConn = true;
else
    isConn = false;
end
