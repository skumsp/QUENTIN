function [] = saveMJNet(Data,outfile)
    Sequences_char = char(Data.Sequence);
    [SeqNum, PosNum] = size(Sequences_char);
    [P_info, Q1, Map] = calculate_P_info(Sequences_char);
    Sequences_MJN = Sequences_char(:, Map);
    Names_MJN = Data.Header;
    delete(outfile);
    fid = fopen(outfile, 'w');
    fprintf(fid, '  ;1.0\n');

    for v = 1:length(Map)%PosNum
      fprintf(fid, '%d;', Map(v));
    end    
    fprintf(fid, '\n');
    for v = 1:length(Map)%PosNum
      fprintf(fid, '10;');
    end    
    fprintf(fid, '\n');
    for e = 1:SeqNum-1
        [token,remain] = strtok(Data(e).Header,'_');
        t = [token '_' int2str(e)];
      fprintf(fid, '>%s', t);
      fprintf(fid, ';1;2;;;;;;\n');
      fprintf(fid, '%s\n', upper(Sequences_MJN(e, :)));
    end
    e = SeqNum;
    [token,remain] = strtok(Data(e).Header,'_');
    t = [token '_' int2str(e)];
    fprintf(fid, '>%s', t);
    fprintf(fid, ';1;2;;;;;;\n');
    fprintf(fid, '%s', upper(Sequences_MJN(e, :)));
    status = fclose(fid)