inputFolder = 'AA_clipped';
splitsTreeFolder = 'SplitsTree';
distThr = 1100;
nIterSimul = 3000;
nIterMCMC = 250;
nInstMCMC = 4;
[HostNetThr, sources, transNets,TransTrees] = quentin(inputFolder,splitsTreeFolder,[],distThr,[],nIterSimul,nIterMCMC,nInstMCMC,[],[]);
sources

bg = biograph(HostNetThr,[],'ShowWeights','on');
view(bg);

TransNet = transNets{1};
bg = biograph(sparse(TransNet),[],'ShowWeights','on');
view(bg)