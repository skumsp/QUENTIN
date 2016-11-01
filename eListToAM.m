function AM = eListToAM(E,nVert)
AM = zeros(nVert, nVert);
nE = size(E,2);
for i=1:nE
    e = E{i};
    if size (e,2) == 2
        AM(e(1),e(2)) = 1;
    else
       AM(e(1),e(2)) = e(3);
    end
end