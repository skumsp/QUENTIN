Welcome to the QUENTIN wiki!

This repository contains the version of QUENTIN implemented in matlab. QUENTIN is also available in python, see

https://github.com/walkergussler/quentin

The description below is about matlab version. The main script is called quentin.m Input/output format:

[HostNetThr, sources, transNets,transTrees] = quentin(inputFolder,splitsTreeFolder,distType,distThr,nclust,nIterSimul,nIterMCMC,nInstMCMC,interHostCoeffs,rho)
Input parameters:

inputFolder is a folder with input files in fasta format. All sequenced should be aligned
splitsTreeFolder is a folder, where SplitsTree is located (could be downloaded at http://www.splitstree.org/)
distType - type of genetic distance used by the algorithm. Possible types: [] (default value), 'evol' (simulation-based distance described in the paper), 'consensus' (distance between consensuses), 'mindist' (minimal distance between populations). Default value: 'evol'
distThr - distance threshold for transmission clusters detection
nclust - number of intra-host variants in each population, to which all populations are reduced to eliminate sampling bias. If nclust = [], then the default value nclust = 12 is used
nIterSimul - maximal number of iterations for viral evolution simulation
nIterMCMC - maximal number of iterations for MCMC algorithm
nInstMCMC - number of parallel instances of MCMC
interHostCoeffs - range of values for the coefficient (\alpha)^{-1} for the estimation of likelihood of genetic distances. If interHostCoeffs = [], then the default value interHostCoeffs = [1 1.25 1.5] is used
rho - range of values for the coefficient \rho for prior probability of transmission tree estimation. If rho = [], then the default value rho = 1 is used
Output:

HostNetThr - graph obtained from the host network by removal of arcs with weights exceeding the threshold. Weakly connected components of this graph are transmission clusters
sources - sources of outbreaks for all transmission clusters
transNets - transmission networks for all transmission clusters
transTrees - transmission trees for all transmission clusters
Example:

inputFolder = 'AA_clipped';

splitsTreeFolder = 'SplitsTree';

distThr = 1100;

nIterSimul = 3000;

nIterMCMC = 250;

nInstMCMC = 4;

[HostNetThr, sources, transNets,TransTrees] = quentin(inputFolder,splitsTreeFolder,[],dis
