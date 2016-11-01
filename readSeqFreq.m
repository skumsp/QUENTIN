function [Sequences_fasta,freq,subtype] = readSeqFreq(seq_file,freq_thr,ntoktype)

           Sequences_fasta = fastaread(seq_file);
           [pathstr,name,ext] = fileparts(seq_file);
           names = char (Sequences_fasta.Header);
           [name, remainGenot] = strtok(name, '_');
          
           nseq = size(names,1);
%            if nseq == 1
%                continue;
%            end
           freq = zeros(1,nseq);
           for u=1:nseq
               s = names(u,:);
               [token,remain] = strtok(s,'_');
               while ~isempty(remain)
                   s = remain;
                   [token,remain] = strtok(s,'_');
               end
               freq(u) = str2num(token);
           end
           
           
           lowFreq = find(freq <= freq_thr);
           if size(lowFreq,2) ~= nseq
               Sequences_fasta(lowFreq) = [];
               freq(lowFreq) = [];
               nseq = size(freq,2);
           end
           
           for j=2:ntoktype
               [next,remainGenot] = strtok(remainGenot, '_.');
           end
           subtype = next;