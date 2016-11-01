function [centroids, freq_centr] = calc_centroids(seq_fasta, freq,nclust, idclustmethod,idtreemethod,treetype, name_pat)

       names = char (seq_fasta.Header);
       nseq = size(seq_fasta,1);
              
       if nseq >= nclust           
           Sequences_fasta_mult = [];
           for u=1:size(seq_fasta,1)
               for v=1:freq(u)
                   seq = seq_fasta(u);
                   seq.Header = [seq.Header, '_', int2str(v)];
                   Sequences_fasta_mult = [Sequences_fasta_mult; seq];
               end
           end
           seqs = Sequences_fasta_mult;
           names = char (seqs.Header);
           nseq = size(Sequences_fasta_mult,1);
           freq = ones(1,nseq);


           DM = seqpdist(seqs,'method','jukes-cantor','indels','pairwise-delete');
           nseq = size(seqs,1);
           
           if strcmp(treetype,'linkage')
                Z = linkage(DM,idtreemethod);
    % %             dendrogram(Z,0);
                T = cluster(Z,'maxclust',nclust);
           end
           if strcmp(treetype,'seqneighjoin')
               names_tree = arrayfun(@num2str, 1:nseq, 'unif', 0);
               Z = seqneighjoin(DM,idtreemethod,names_tree);                 

               [LeafClusters, NodeClusters] = cluster(Z,[],'criterion',idclustmethod,'maxclust',nclust);

    %            h = plot(Z,'orient','top');
    %            set(h.BranchLines(NodeClusters==4),'Color',[0 0 0])
    %            set(h.BranchLines(NodeClusters==3),'Color',[0 1 0])
    %            set(h.BranchLines(NodeClusters==2),'Color',[0 0 1])
    %            set(h.BranchLines(NodeClusters==1),'Color',[1 0 0])
    %            for jj=5:nclust
    %                 set(h.BranchLines(NodeClusters==jj),'Color',[rand rand rand])
    %            end

               nclust = max(LeafClusters);
               leafs = get(Z,'LeafNames');
               T = zeros(nseq,1);
               for u = 1:nseq
                    T(str2num(leafs{u})) = LeafClusters(u);
               end
           end

           centroids = [];
           freq_centr = [];
          
           hm = containers.Map;
           for c=1:nclust
               clust = find(T == c);
               seqs_clust = [];
               for u=1:size(clust,1)
                   for v=1:freq(clust(u))
                       seqs_clust = [seqs_clust; seqs(clust(u))];
                   end
               end
               cons = seqconsensus(seqs_clust);
               f = sum(freq(clust),2);
               if isKey(hm,cons)
                   fprev = hm(cons);
                   hm(cons) = fprev + f;
               else
                   hm(cons) = f;
               end
           end
           cens = keys(hm);
           for c=1:size(cens,2)		
                nucleo = cens{c}; 
                f = hm(nucleo);
                centr.Header = [name_pat,'_cen_',int2str(f),'_',int2str(c)];
                centr.Sequence = nucleo;
                freq_centr = [freq_centr; f];
                centroids = [centroids; centr];
           end
       else
           centroids = [];
           freq_centr = [];
           for u=1:nseq
                centr.Header = [name_pat,'_cen_',int2str(freq(u)),'_',int2str(u)];
                centr.Sequence = seq_fasta(u).Sequence;
                freq_centr = [freq_centr; freq(u)];
                centroids = [centroids; centr];
           end
       end