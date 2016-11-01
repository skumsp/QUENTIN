function alignment = align_IgnoreGaps(seqs)
    if size(seqs,1) >= 3
        alignment = multialign(seqs);
%         alignment = multialign(seqs,'GapOpen', 15,'ExtendGap',6);
    else
        alignment = seqs;
        return;
    end
    seqchar = char (alignment.Sequence);
    l = size(seqchar,2);
    nseq = size(seqchar,1);
    gapsmember = zeros(1,l);
    for i=1:l
        gapsmember(i) = ismember('-',seqchar(:,i));
    end
    gaps = find(gapsmember == 1);
    seqchar(:,gaps) = [];
    for i=1:nseq
        alignment(i).Sequence = seqchar(i,:);
    end
