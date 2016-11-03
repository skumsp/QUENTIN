function [d12, d21] = outbreaks_aspen1(seq1,freq1,subtype1,name1,p1,seq2,freq2, subtype2,name2,p2,outbName,distType,timeInter,splitsTreeFolder)

%     timeInter = 5000;
    epsilon = 0.0001
    
    toAlign = true;
    toIgnoreGaps = true;

   
    
    if strcmp(name1,name2) == 1
        d12 = 1;
        d21 = 1;
        resDir = [outbName '_results'];
%         resDir = [outbName '_results_' distType];
        outfile = [resDir filesep int2str(p1) '_' int2str(p2) '_' name1 '_' name2 '_res.txt'];
        fid = fopen(outfile, 'w');
        fprintf(fid, [num2str(d12) '\n']);
        fprintf(fid, [num2str(d21) '\n']);
        status = fclose(fid);
        return;
    end
    
%     if strcmp(subtype1,subtype2) == 0
%         d12 = intmax;
%         d21 = intmax;
%         resDir = [outbName '_results'];
% %         resDir = [outbName '_results_' distType];
%         outfile = [resDir filesep int2str(p1) '_' int2str(p2) '_' name1 '_' name2 '_res.txt'];
%         fid = fopen(outfile, 'w');
%         fprintf(fid, [num2str(d12) '\n']);
%         fprintf(fid, [num2str(d21) '\n']);
%         status = fclose(fid);
%         return;
%     end
    
    if strcmp(distType,'mindist')
        seq1 = char (seq1.Sequence);
        seq2 = char (seq2.Sequence);
        dm = pdist2(seq1, seq2,'hamming');
        mindist = min(min(dm));
        if mindist == 0
            mindist = epsilon;
        end
%         for i = 1:size(seq1,1)
%             for j = 1:size(seq2,1)
%                 d = pdist2(seq1(i,:), seq2(i,:),'hamming');
%                 if d < mindist
%                     mindist = d;
%                 end
%             end
%         end
        d12 = mindist;
        d21 = mindist;
        resDir = [outbName '_results'];
        outfile = [resDir filesep int2str(p1) '_' int2str(p2) '_' name1 '_' name2 '_res.txt'];
        fid = fopen(outfile, 'w');
        fprintf(fid, [num2str(d12) '\n']);
        fprintf(fid, [num2str(d21) '\n']);
        status = fclose(fid);
        return;
    end
    
    if strcmp(distType,'consensus')
        cons1 = seqconsensus(seq1);
        cons2 = seqconsensus(seq2);
        dm = pdist2(cons1, cons2,'hamming');
        d12 = dm;
        d21 = dm;
        resDir = [outbName '_results_' distType];
        outfile = [resDir filesep int2str(p1) '_' int2str(p2) '_' name1 '_' name2 '_res.txt'];
        fid = fopen(outfile, 'w');
        fprintf(fid, [num2str(d12) '\n']);
        fprintf(fid, [num2str(d21) '\n']);
        status = fclose(fid);
        return;
    end

    
%     Centroids_fasta_all = [seq1; seq2];    
%     delete(strcat(outbName,'_centr.fasta'))
%     fastawrite(strcat(outbName,'_centr.fasta'),Centroids_fasta_all);
     
%     if strcmp(subtypes_pat{1},subtypes_pat{2})
%       [d12, d21] = calcTimeSamp_MJNSimul(seq1,seq2,freq1,freq2,name1,name2,splitsTreeFolder, timeInter, toAlign, toIgnoreGaps);
      try
        [d12, d21] = calcTimeSamp_MJNSimul(seq1,seq2,freq1,freq2,name1,name2,splitsTreeFolder, timeInter, toAlign, toIgnoreGaps);
      catch 
          name1
          name2
          error('some error')         
      end
      
    fileMJ_in = ['MJ_' name1 '_' name2 '_in.nex'];
    fileMJ_out = ['MJ_' name1 '_' name2 '_out.nex'];
    delete(fileMJ_in);
    delete(fileMJ_out);
