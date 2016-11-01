function pids = parsePIDsWin(pidsAr, proc)
    indCom =  find(strcmp(proc, pidsAr));
    if size(indCom,2) == 0
        pids = {};
        return;
    end
    indPID = indCom + 1;
    pids = pidsAr(indPID);