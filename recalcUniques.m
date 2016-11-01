function [seqNew,freqNew] = recalcUniques(seqs,freqs)
           seqNew = [];
           freqNew = [];
           hm = containers.Map;
           fr = containers.Map;
           for i=1:size(seqs,1)
               nucl = seqs(i).Sequence;
               name = seqs(i).Header;
               f = freqs(i);

               if ~isKey(hm,nucl)
                   hm(nucl) = name;
                   fr(nucl) = f;
               else
                   f_old = fr(nucl);
                   fr(nucl) = f_old + f;
               end
           end
           nucls = keys(hm);
           for c=1:size(nucls,2)		
                nucleo = nucls{c}; 
                name = hm(nucleo);
                f = fr(nucleo);
                centr.Header = name;
                centr.Sequence = nucleo;
                seqNew = [seqNew; centr];
                freqNew = [freqNew; f];
           end