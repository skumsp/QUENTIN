function [HostNetThr, sources, transNets,transTrees] = quentin(inputFolder,splitsTreeFolder,distType,distThr,nclust,nIterSimul,nIterMCMC,nInstMCMC,interHostCoeffs,rho)

%todo 1) return networks for every component 2) linkage vs neighbor joining

outdir = inputFolder;
if isempty(distType)
    distType = 'evol';
end
% distType = 'consensus';
% distType = 'mindist';
ntoktype = 2;
freq_thr = 0;
if isempty(nclust)
    nclust = 12;
end

files = dir(fullfile(outdir,'*.fas'));
nSamp = size(files,1);
names_pat = cell(1,nSamp);
centroids_pat = cell(1,nSamp);
centroids_pat_freq = cell(1,nSamp);

for i =1:size(files,1) 
   [pathstr,name,ext] = fileparts(files(i).name);
   namelong = name;
   ['Calculating clusters ' namelong]   
      name = strtok(name, '_');
      seq_file = [outdir,filesep,files(i).name];   
      [Sequences_fasta,freq,subtype] = readSeqFreq(seq_file,freq_thr,ntoktype);
    if strcmp(distType,'evol') 
      [centroids, freq] = calc_centroids(Sequences_fasta, freq, nclust, [],'weighted','linkage',name);
    else
        centroids = Sequences_fasta;
    end   
   names_pat{i} = name;
   centroids_pat{i} = centroids;
   centroids_pat_freq{i} = freq;
   
end

[X,Y] = ndgrid(1:nSamp,1:nSamp);
product = [X(:) Y(:)];

DSampPar = zeros(1,size(product,1));
DSampParRev = zeros(1,size(product,1));

x = product(:,1);
y = product(:,2);

for p=1:size(product,1)
        p1 = x(p);
        p2 = y(p);
        if p1>=p2
            continue;
        end

        try
            [d12, d21] = outbreaks_aspen1(centroids_pat{p1},centroids_pat_freq{p1},[],names_pat{p1},p1,centroids_pat{p2},centroids_pat_freq{p2},[],names_pat{p2},p2,inputFolder,distType,nIterSimul,splitsTreeFolder);
        catch
            fid = fopen( 'log.txt', 'w' );
            message = [names_pat{p1},' ', names_pat{p2}];
            fprintf(fid, '%s\n', message);	
            fclose(fid);
            exit(1);
        end
        DSampPar(p) = d12;
        DSampParRev(p) = d21;
end

for i=1:size(product,1)
        p1 = product(i,1);
        p2 = product(i,2);
        if p1>=p2
            continue;
        end
        DSamp(p1,p2) = DSampPar(i);
        DSamp(p2,p1) = DSampParRev(i);
end


AMSamp_dir = (DSamp <= DSamp');
AMSamp = AMSamp_dir.*(DSamp <=distThr);
HostNetThr = DSamp.*AMSamp;

[V,D] = eig(AMSamp);
eigval = diag(D);
k = find(eigval == max(eigval));
centr = V(:,k)/sum(V(:,k),1);
npat = size(AMSamp,2);

res = zeros(3,npat);
res(1,:) = 1:npat;
res(2,:) = centr;
res = sortrows(res',-2)'
names_pat(res(1,:))

source = find(centr == max(centr));
if size(source,1) == 2
    if DSamp(source(1),source(2)) < DSamp(source(2),source(1))
        source = source(1);
    else
        if DSamp(source(2),source(1)) < DSamp(source(1),source(2))
            source = source(2);
        end
    end
end

sources = names_pat(source);

[S, C] = graphconncomp(sparse(AMSamp),'Weak',true);

maxdist = [];
nEdgeModif = 2;
if isempty(rho)
    rho = [1];
end
ndecr = 4;
if isempty(interHostCoeffs)
    interHostCoeffs = [1 1.25 1.5];
end

transNets= cell(1,S);
transTrees = cell(1,S);
for c=1:S
   comp = find(C == c);
   nComp = size(comp,2);
   if nComp > 1
       DSamp_comp = DSamp(comp,comp);
       
%        [transNetsComp,TransTreesComp,record,aux] = findTransNetMCMC4(DSamp_comp,nIterMCMC,maxdist,nEdgeModif,rho,ndecr,[], 0,[],[]);
       
       transNetsMCMCPar = cell(1,nInstMCMC);
       TransNetTreesWeightPar = cell(1,nInstMCMC);
       recordPar = zeros(1,nInstMCMC);
       parfor i=1:nInstMCMC
            interHostCoeff = datasample(interHostCoeffs,1);
            [transNetsMCMCPar{i},TransNetTreesWeightPar{i},recordPar(i),aux] = findTransNetMCMC5(DSamp_comp,nIterMCMC,maxdist,nEdgeModif,rho,interHostCoeff,ndecr,[], 0,[],[]);
       end
       [argvalue, argmax] = max(recordPar);
       transNetsComp = transNetsMCMCPar{argmax};
       TransTreesComp = TransNetTreesWeightPar{argmax};
       
       transNets{c} = transNetsComp{1};
       transTrees{c} = TransTreesComp{1};
   end
end

['finshed']
