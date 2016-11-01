function [d12, d21] = calcTimeSamp_MJNSimul(seqPat1,seqPat2,freq1,freq2,namePat1,namePat2,splitsTreeFolder,timeInter, toAlign, toIgnoreGaps)
    warning('error', 'MATLAB:nearlySingularMatrix');
    mutprob = 0.01;

    cutoff_heter = 0;
    epsilonMJN = 0;
    maxPopl = 10^(12);
    maxTimeRun = 30;
    maxAttempts = 2;
    tpause = 1;

            nseq1 = size(seqPat1,1);
            nseq2 = size(seqPat2,1);

            seqBoth = [seqPat1; seqPat2];
            if toAlign
                seqBoth = align(seqBoth,toIgnoreGaps);
                seqPat1_align = seqBoth(1:nseq1);
                seqPat2_align = seqBoth((nseq1+1):end);
                [seqPat1_new,freq1_new] = recalcUniques(seqPat1_align,freq1);
                [seqPat2_new,freq2_new] = recalcUniques(seqPat2_align,freq2);
                seqBoth = [seqPat1_new; seqPat2_new];
                nseq1 = size(seqPat1_new,1);
                nseq2 = size(seqPat2_new,1); 
                freq1 = freq1_new;
                freq2 = freq2_new;
            end
            
            
%             delete('align.fasta');
%             fastawrite('align.fasta',seqBoth);
            
            [P_info, Q1, Map] = calculate_P_info(char(seqBoth.Sequence),cutoff_heter);
            len_eff = size(Map,2);
            fileMJ_in = ['MJ_' namePat1 '_' namePat2 '_in.nex'];
            fileMJ_out = ['MJ_' namePat1 '_' namePat2 '_out.nex'];
            
            
            tryST = true;
            attempts = 1;
            while tryST
                tryST = false;
                delete(fileMJ_in);
                delete(fileMJ_out);
                saveNEXUSSplitsTree(seqBoth,fileMJ_in,cutoff_heter, epsilonMJN);
                
%                 Linux pid
%                 compid = 'pidof java';
%               Windows pid
                compid = 'tasklist /fi "Imagename eq SplitsTree.exe"';
                [status, pids] = system(compid);
                pidsOutput = strsplit(pids);
                pidsAr = parsePIDsWin(pidsOutput,'SplitsTree.exe');
                
%                 newsplitsTreeFolder = [splitsTreeFolder namePat1 '_' namePat2];
%                 copyfile(splitsTreeFolder,newsplitsTreeFolder);

                command = [splitsTreeFolder filesep 'SplitsTree +g -x ' char(39) 'Execute FILE=' fileMJ_in '; Update; SAVE FILE=' fileMJ_out '; QUIT' char(39) '&'];      
    %             status = system(command);
                ['Running SplitsTree']
                status = system(command);

                pause on
                pause(1);
                [status, pidsnew] = system(compid);
                pidsOutputNew = strsplit(pidsnew); 
                pidsArNew = parsePIDsWin(pidsOutputNew,'SplitsTree.exe');
                %Windows
                currPID = setdiff(pidsArNew,pidsAr);
                %Linux
%                 currPID = setdiff(pidsOutputNew,pidsOutput);
                pause off;
                
                pause on
                iwait = 1;
                while exist(fileMJ_out,'file') == 0;
                    ['Waiting for SplitsTree: ' int2str(iwait)]
                    pause(1);
                    iwait = iwait + 1;
                    if iwait > maxTimeRun
                        ['Time exceeded']
                        currPID = currPID{1}
%                         Linux
%                         comkill = ['kill ' currPID];
%                       Windows
                        comkill = ['taskkill /f /pid ' currPID];
                        system(comkill);
                        tryST = true;
                        attempts = attempts + 1;
                        epsilonMJN = epsilonMJN + 1;
                        if attempts > maxAttempts
                            d12 = -1;
                            d21 = -1;
%                             rmdir(newsplitsTreeFolder,'s');
                            return;
                        else
                            break;
                        end
                        
                    end
                end
                pause off
            end
            
            
            if status == 0
                ['MJN calculated']
                pause(tpause);
            end

    %         fastawrite(fileMJ,seqBoth);
            try
                [D, eqVert] = readNEXUS2(fileMJ_out);
            catch 
                warning('SplitsTree crashed');
                d12 = -1;
                d21 = -1;
%                 rmdir(newsplitsTreeFolder,'s');
                return;
            end
            D_all = D;
            nseq2 = nseq2 - eqVert;

            SM_all = zeros(size(D_all,1),size(D_all,2));
            for u=1:size(D_all,1)
                for v=(u+1):size(D_all,2)
                    SM_all(u,v) = D_all(u,v)*(mutprob/3)^(D_all(u,v))*(1-mutprob)^(len_eff - D_all(u,v));
                    SM_all(v,u) = SM_all(u,v);
                end
                SM_all(u,u) = (mutprob/3)^(D_all(u,u))*(1-mutprob)^(len_eff - D_all(u,u));
            end

            Q_all = SM_all;
            nseq_all = size(Q_all,2);

%             covtime = simulEvol(Q_all,nseq1,nseq2,timeInter,maxPopl);
            covtime = simulEvol1(Q_all,nseq1,nseq2,timeInter,maxPopl);
%             immRespStr = 0.8;
%             immRespRate = 0.04;
%             immRespDecay = 0.04;
%             covtime = simulEvolImmResp(Q_all,freq1,freq2,nseq1,nseq2,timeInter); 
            d12 = covtime;

            perm = [nseq1+1:nseq1+nseq2 1:nseq1 nseq1+nseq2+1:nseq_all];
            Q_all = Q_all(perm,perm);  
            nseq2 = nseq2 + eqVert;
            nseq1 = nseq1 - eqVert;
%             covtime = simulEvol(Q_all,nseq2,nseq1,timeInter,maxPopl); 
            covtime = simulEvol1(Q_all,nseq2,nseq1,timeInter,maxPopl); 
%             covtime = simulEvolImmResp(Q_all,freq2,freq1,nseq2,nseq1,timeInter); 
            d21 = covtime;
%             rmdir(newsplitsTreeFolder,'s');



