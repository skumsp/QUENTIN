clear;
splitsTreeFolder = 'SplitsTree';
outbName = 'Vietnam_1a';
% outdir = [outbName,' ','outbreak'];
outdir = outbName;
distType = 'evol';
% distType = 'consensus';
% distType = 'mindist';
ntoktype = 2;
freq_thr = 0;
nclust = 12;
timeInter = 100000;

thr = 1100;
% thr = 0.06;
% thr = 0.0377;

resDir = [outbName '_results'];
mkdir(resDir);
centDir = [outbName '_centroids'];
mkdir(centDir);


files = dir(fullfile(outdir,'*.fas'));
nSamp = size(files,1);
names_pat = cell(1,nSamp);
subtypes_pat = cell(1,nSamp);
centroids_pat = cell(1,nSamp);
centroids_pat_freq = cell(1,nSamp);
sequences_fasta_all = [];

parfor i =1:size(files,1) 
   [pathstr,name,ext] = fileparts(files(i).name);
   namelong = name;
   ['Calculate clusters ' namelong]
   
   outfile_cent = [centDir filesep namelong '_cent.mat'];
   if exist(outfile_cent,'file') == 2
       load(outfile_cent);
   else
      name = strtok(name, '_');
      seq_file = [outdir,filesep,files(i).name];   
      [Sequences_fasta,freq,subtype] = readSeqFreq(seq_file,freq_thr,ntoktype);
%       sequences_fasta_all = [sequences_fasta_all; Sequences_fasta];
    if strcmp(distType,'evol') 
      [centroids, freq] = calc_centroids(Sequences_fasta, freq, nclust, [],'weighted','linkage',name);
%       m=matfile(outfile_cent,'writable',true)
%       m.name=name;
%       m.centroids=centroids;
%       m.subtype = subtype;
%       m.freq = freq;
    else
        centroids = Sequences_fasta;
    end
   end
   
   names_pat{i} = name;
   subtypes_pat{i} = subtype;
   centroids_pat{i} = centroids;
   centroids_pat_freq{i} = freq;
   
%    save(outfile_cent,'name','subtype','centroids','freq_centr');
end

% fastawrite([outbName '_all.fas'],sequences_fasta_all);

[X,Y] = ndgrid(1:nSamp,1:nSamp);
product = [X(:) Y(:)];

DSampPar = zeros(1,size(product,1));
DSampParRev = zeros(1,size(product,1));

x = product(:,1);
y = product(:,2);

if exist(resDir,'file') == 0
     mkdir(resDir);
end

parfor p=1:size(product,1)
        p1 = x(p);
        p2 = y(p);
        if p1>=p2
            continue;
        end
        

        outfile = [resDir filesep int2str(p1) '_' int2str(p2) '_' names_pat{p1} '_' names_pat{p2} '_res.txt'];
        if exist(outfile,'file') == 2
            continue;
            M = csvread(outfile);
            d1 = M(1);
            if d1 ~= -1
                continue;
            end
        end
        
%          [d12, d21] = outbreaks_aspen(int2str(p1),int2str(p2),outbName,distType,splitsTreeFolder);
        try
            [d12, d21] = outbreaks_aspen1(centroids_pat{p1},centroids_pat_freq{p1},subtypes_pat{p1},names_pat{p1},p1,centroids_pat{p2},centroids_pat_freq{p2},subtypes_pat{p2},names_pat{p2},p2,outbName,distType,timeInter,splitsTreeFolder);
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

% timeInter = 2500;
% delta = 0;
% AMSamp = zeros(nSamp,nSamp);
% 
% for u=1:nSamp
%     for v=(u+1):nSamp
%             m = min(DSamp(u,v),DSamp(v,u));
%             if m<=thr
%                 if DSamp(u,v) == m
%                     AMSamp(u,v) = 1;
%                 end
%                  if DSamp(v,u) == m
%                     AMSamp(v,u) = 1;
%                 end
%             end
%     end
% end

AMSamp_dir = (DSamp <= DSamp');
AMSamp = AMSamp_dir.*(DSamp <=thr);
DSamp_dir = DSamp.*AMSamp;

bg = biograph(DSamp_dir,names_pat,'ShowWeights','on');
view(bg);

% U = DSamp < 3000;
% bg = biograph(DSamp.*U,names_pat,'ShowWeights','on');
% view(bg);

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

sources = names_pat(source)

% nshow = 4;
% if nshow < nSamp
%     toshow = res(1,1:nshow);
%     bg = biograph(DSamp_dir(toshow,toshow),names_pat(toshow),'ShowWeights','on');
%     view(bg);
% end

%% 

[S, C] = graphconncomp(sparse(AMSamp),'Weak',true);
colors = {};
realSource = 'AI004';
for i=1:nSamp
    if strcmp(names_pat{i},realSource)
        colors{i} = 'Red';
    else
        colors{i} = 'Blue';
    end
end

% nonIsol= [];
% for i=1:nSamp
%     comp = find(C == C(i));
%     if size(comp,2) > 1
%         nonIsol = [nonIsol i];
%     end
% end

% AMSampPaj = AMSamp_dir(nonIsol,nonIsol);
% names_patPaj = names_pat(nonIsol);
% colorsPaj = colors(nonIsol);

adj2pajek(AMSamp,names_pat,colors,[],[outbName '_' distType '_pajek.net'])
%% 

% transNet = findTransNetSPT(DSamp_dir,0);
% bgTrans = biograph(transNet,names_pat,'ShowWeights','on');
% view(bgTrans);

% transNetSPT = findTransNetSPT(DSamp_dir,1);
% bgTrans = biograph(transNetSPT,names_pat,'ShowWeights','on');
% view(bgTrans);
% wSPT = sum(sum(transNetSPT,2),1)
% %% 
% transNetArb = findTransNetArb(DSamp_dir,'min');
% bgTrans = biograph(transNetArb,names_pat,'ShowWeights','on');
% view(bgTrans);
% wArb = sum(sum(transNetArb,2),1)
% %% 
% transNetQIP = findTransNetQIP(DSamp_dir,source);
% bgTrans = biograph(transNetQIP,names_pat,'ShowWeights','on');
% view(bgTrans);
%% 
nIter = 250;
debugMode = 0;
maxdist = 3000;
nEdgeModif = 2;
% pow = 1;
powers = [1];
pow = 1;
ndecr = 4;

for c=1:S
   comp = find(C == c);
   nComp = size(comp,2);
   if nComp > 1
%        AMSamp_comp = AMSamp(comp,comp);
%        bg = biograph(AMSamp_comp,names_pat(comp),'ShowArrows', 'on');
%        view(bg);
       
       DSamp_comp = DSamp(comp,comp);
% add type of init tree
         [transNetsMCMC,TransNetTreesWeight,record,aux] = findTransNetMCMC4(DSamp_comp,nIter,maxdist,nEdgeModif,powers,ndecr,[], 0,[],[]);
%        [transNetsMCMC,TransNetTreesWeight,record] = findTransNetMCMC(DSamp_comp,nIter,debugMode);
%        [transNetsMCMC,TransNetTreesWeight,record,aux] = findTransNetMCMC3(DSamp_comp,nIter,maxdist,nEdgeModif,pow,[], 0,[],[]);
%        [transNetsMCMC,TransNetTreesWeight,record,auxc] = findTransNetMCMC3(DSamp,nIter,maxdist,nEdgeModif,pow,[],0,[],aux);

       for i=1:size(transNetsMCMC,2)
            TransNet = transNetsMCMC{i};
            bg = biograph(sparse(TransNet),names_pat(comp),'ShowWeights','on');
            view(bg);
       end
%        heights = zeros(1,size(TransNetTreesWeight,2));
%        for i=1:size(TransNetTreesWeight,2)
%             TransNetTree = TransNetTreesWeight{i};
%             [DM_tree,labels] = treeDMdirect(TransNetTree,names_pat(comp));
%             bg = biograph(sparse(DM_tree),labels,'ShowWeights','on');
%             view(bg);
%             heights(i) = getHeight(TransNetTree);
%        end
   end
end

['finshed']
%% 

% indeg = sum(AMSamp,1);
% K = diag(indeg) - AMSamp;
% K(source,:) = [];
% K(:,source) = [];
% det(K)
% 
% gamma = 2.5;
% bound = 1500;
% transNetBaB = findTransNetBaB(DSamp_dir,gamma,source,bound);
% bgTrans = biograph(transNetBaB,names_pat,'ShowWeights','on');
% view(bgTrans);
%% 

% timepoints = [0 1.8 2.8 3.3 4.6 6.0 7.2 7.6];
% nTimePoints = 8;
% res = zeros((nTimePoints)*(nTimePoints-1)/2,2);
% pairs = 1;
% add = 0;
% for i = 1:nTimePoints
%     for j = (i+1):nTimePoints
%         IndexC = strfind(names_pat, ['RL',int2str(i)])
%         ii = find(not(cellfun('isempty', IndexC)))
%         IndexC = strfind(names_pat, ['RL',int2str(j)])
%         jj = find(not(cellfun('isempty', IndexC)))
%         res(pairs,1) = abs(timepoints(j)-timepoints(i));
%         res(pairs,2) = min(DSamp(ii,jj),DSamp(jj,ii));
%         pairs = pairs + 1;
%     end
% end
