function AM = eListToAM_matr(E,nVert,symm)
AM = zeros(nVert, nVert);
nE = size(E,1);
for i=1:nE
    e = E(i,:);
    if size (e,2) == 2
        AM(e(1),e(2)) = 1;
        if symm
           AM(e(2),e(1)) = 1; 
        end
    else
       AM(e(1),e(2)) = e(3);
       if symm
           AM(e(2),e(1)) = e(3); 
       end
    end
end