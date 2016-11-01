clear;
outbName = 'all_clipped';
% outbName = 'IN_clean_all_Candidate_singletons_removed';
outdir = outbName;
thr = 3000;
ntokgenot = 5;
genAn = '1a';

files = dir(fullfile(outdir,'*.fas'));
nSamp = size(files,1);

isComp = zeros(1,nSamp);
timeClust = zeros(1,nSamp);
mindistClust = zeros(1,nSamp);
names_pat = cell(1,nSamp);

cons = [];

for i =1:nSamp
       [pathstr,name,ext] = fileparts(files(i).name);
       [name,res] = strtok(name, '_');
       name
       for j = 2:ntokgenot
           [gen,res] = strtok(res, '_');
       end
%        if ~strcmp(gen,genAn)
%            continue;
%        end
       names_pat{i} = files(i).name;
       seq_file = [outdir,'\',files(i).name];
       seq = fastaread(seq_file);
       CSeq = seqconsensus(seq);
       con.Header = name;
       con.Sequence = CSeq;
       cons = [cons; con];
end
% outfile = [outdir,'\','consIndiana','_',genAn,'.fas'];
outfile = ['cons_',outbName,'.fas'];
fastawrite(outfile,cons);