function v = saveparfor(filename,res)
%     filename = ['res_',int2str(t),'.mat'];
    save(filename,'res');
    v = 1;