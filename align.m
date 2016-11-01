function alignment = align(seqs,toIgnoreGaps)
    if size(seqs,1) >= 3
%         alignment = multialign(seqs);
        alignment = multialign(seqs,'GapOpen', 15,'ExtendGap',6);
    else
        [Score, alignment] = nwalign(seqs(1),seqs(2),'GapOpen', 15,'ExtendGap',6);
        alSeq1.Header = seqs(1).Header;
        alSeq2.Header = seqs(2).Header;
        alSeq1.Sequence = alignment(1,:);
        alSeq2.Sequence = alignment(3,:);
        alignment = [alSeq1; alSeq2];
    end
    seqchar = char (alignment.Sequence);    
    l = size(seqchar,2);
    nseq = size(seqchar,1);
    gapsmember = zeros(1,l);
    for i=1:l
        gapsmember(i) = ismember('-',seqchar(:,i));
    end
    gaps = find(gapsmember == 1);
    if toIgnoreGaps
        seqchar(:,gaps) = [];
    else
        stpos = 1;
        endpos = l;
        while gapsmember(stpos) == 1
            stpos = stpos+1;
        end
        while gapsmember(endpos) == 1
            endpos = endpos-1;
        end
        seqchar = seqchar(:,stpos:endpos);
    end
    for i=1:nseq
        alignment(i).Sequence = seqchar(i,:);
    end
