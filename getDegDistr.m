function distrib = getDegDistr(DM)
AM = DM > 0;
nSamp = size(AM,2);
outdeg = sum(AM,2);
bins = 0:nSamp;
[distrib,bounds] = histcounts(outdeg,bins,'Normalization', 'probability');