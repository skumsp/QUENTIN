function [D,eqVert] = readNEXUS(addr, posSource)
    fid = fopen(addr);
    line = textscan(fid,'%s',1);
    while ~strcmp(line{1},'DIMENSIONS')
        line = textscan(fid,'%s',1);
    end
    line = textscan(fid, '%s%u',1, 'Delimiter', '=;');
    nseq = line{2};
    line = textscan(fid,'%s',1);
    namesSeq = cell(1,nseq);
    for i=1:nseq
        line = textscan(fid,'%s%s',1,'Delimiter', strcat(char(39)));
        namesSeq{i} = line{2};
    end
    line = textscan(fid,'%s',1);
    while ~strcmp(line{1},'Network;')
        line = textscan(fid,'%s',1);
    end
    textscan(fid,'%s',1);
    line = textscan(fid,'%s%u',3,'Delimiter','=');
    mapSeqToNet = zeros(1,nseq);
    nVert = line{2}(2);
    nEdges = line{2}(3);
    line = textscan(fid,'%s',1);
    while ~strcmp(line{1},'TRANSLATE')
        line = textscan(fid,'%s',1);
    end
    line = fgets(fid);
    while line(1)~=';'
        line = fgets(fid);
        C = strsplit(line);
        id = str2num(C{1});
        for i=2:(size(C,2)-1)
%             name = strrep(C{i},char(39),' ');
            name = C{i};
            if name(end) == ','
                name = name(2:end-2);
            else
                name = name(2:end-1);
            end
            for j=1:nseq
                if strcmp(char(namesSeq{j}),name)
                    mapSeqToNet(j) = id;
                end
            end
        end
    end
%     line = textscan(fid,'%u%s%s',nseq,'Delimiter',strcat(char(39),','));
%     for i=1:nseq
%         id = line{1}(i);
%         name = line{2}{i};
%         for j=1:nseq
%             if strcmp(namesSeq{j},name)
%                 mapSeqToNet(j) = id;
%             end
%         end
%     end
    line = textscan(fid,'%s',1);
    while ~strcmp(line{1},'EDGES')
        line = textscan(fid,'%s',1);
    end
    line = textscan(fid,'%u%u%u%s%f',nEdges,'Delimiter',',=');
    edgeList = zeros(nEdges,3);
    for i=1:nEdges
        u = line{2}(i);
        v = line{3}(i);
        d = 1;
        if size(line{4}{i},1) > 0
            d = line{5}(i);
        end
        edgeList(i,1) = u;
        edgeList(i,2) = v;
        edgeList(i,3) = d;
    end
    edgeList1 = [];
    for i = 1:nEdges
        u = edgeList(i,1);
        v = edgeList(i,2);
        if edgeList(i,3) == 1
            edgeList1 = [edgeList1;  u v];
        else
            nNewVert = edgeList(i,3)-1;
            edgeList1 = [edgeList1; u, nVert+1];
            for j = 1:(nNewVert-1)
                edgeList1 = [edgeList1; nVert+j nVert+j+1];
            end
            edgeList1 = [edgeList1; nVert+nNewVert v];
            nVert = nVert + nNewVert;
        end
    end
    D = zeros(nVert,nVert);
    for i = 1:size(edgeList1,1)
        D(edgeList1(i,1),edgeList1(i,2)) = 1;
        D(edgeList1(i,2),edgeList1(i,1)) = 1;
    end
    
    additVertNum = nVert;
    eqVert = 0;
    for i=1:nVert
        maptoi = find(mapSeqToNet == i);
        if size(maptoi,2) == 2
           j = max(maptoi);
           mapSeqToNet(j) = [];
           eqVert = eqVert+1;
        end
      
    end
    
    mapNonseqToNet = setdiff(1:nVert,mapSeqToNet);
    order = [mapSeqToNet,mapNonseqToNet];
    D = D(order,order);
    fclose(fid);