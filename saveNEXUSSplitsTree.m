function [] = seveNEXUS(Data,outfile,cutoff, epsilon)
    Sequences_char = char(Data.Sequence);
    [SeqNum, PosNum] = size(Sequences_char);
%     [P_info, Q1, Map] = calculate_P_info(Sequences_char,cutoff);
%     Sequences_MJN = Sequences_char(:, Map);
    Sequences_MJN = Sequences_char;
    delete(outfile);
    fid = fopen(outfile, 'w');
    fprintf(fid, '#nexus\n');
    fprintf(fid, '\n');
    fprintf(fid, 'BEGIN Taxa;\n');
    fprintf(fid, ['DIMENSIONS ntax=' int2str(SeqNum) ';\n']);
    fprintf(fid, 'TAXLABELS\n');
    
    for i=1:SeqNum
        fprintf(fid, ['[' int2str(i) '] ' char(39) Data(i).Header char(39) '\n']);
    end
    fprintf(fid,';\n');
    fprintf(fid,'END; [Taxa]\n');
    fprintf(fid,'\n');
    
    fprintf(fid,'BEGIN Characters;\n');
    fprintf(fid,['DIMENSIONS nchar=' int2str(size(Sequences_MJN,2)) ';\n']);
    fprintf(fid,'FORMAT\n');
    fprintf(fid,' datatype=DNA missing=? gap=- symbols="atgc" labels=no transpose=no interleave=no;\n');
    fprintf(fid,'MATRIX\n');
    for i = 1:SeqNum
        fprintf(fid,[lower(Sequences_MJN(i,:)) '\n']);
    end
    fprintf(fid,';\n');
    fprintf(fid,'END; [Characters]\n');
    fprintf(fid,'\n');
    fprintf(fid,'BEGIN st_Assumptions;\n');
    fprintf(fid,'uptodate;\n');
    fprintf(fid,['chartransform=MedianJoining Epsilon = ' int2str(epsilon) ' SpringEmbedderIterations = 2000 LabelEdges = false ShowHaplotypes = false SubdivideEdges = false ScaleNodesByTaxa = false;\n']);
    fprintf(fid,' exclude  no missing;\n');
    fprintf(fid,'autolayoutnodelabels;\n');
    fprintf(fid,'END; [st_Assumptions]\n');
    
    status = fclose(fid);