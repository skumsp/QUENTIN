function time = simulEvol1(Q,nseq1,nseq2,timeInter,maxPopl);

        nseq = size(Q,2);
        thr = 1;

        
%         [Evec,Eval] = eig(Q);
%         Eval = diag(Eval);
        x = zeros(nseq,1);
        x(1:nseq1) = ones(nseq1,1);
        E = eye(nseq,nseq);
        
        Q = gpuArray(Q);
        x = gpuArray(x);
        E = gpuArray(E);
        for t=2:timeInter
            t
            x = (1-sum(x,1)/maxPopl) * (E+Q)*x;
%             x(:,t) = x(:,t-1) + (1-sum(x(:,t-1),1)/maxPopl) *Q*x(:,t-1);
            popl = find(x((nseq1+1):nseq1+nseq2) >= thr);
            if size(popl,1) == nseq2
                break;
            end
        end
        if t < timeInter
            time = t;
        else
            time = intmax;
        end
        
%         totalPopl = sum(x,1);
%         figure
%         plot(totalPopl(1:t));
%         A = (Q > 0) - eye(nseq,nseq);
%         pajekFile = 'MJNpajek.net';
%         colors = cell(1,nseq);
%         names = cell(1,nseq);
%         for i = 1:nseq1
%             colors{i} = 'Red';
%             names{i} = int2str(i);
%         end
%         for i = (nseq1+1):(nseq1+nseq2)
%             colors{i} = 'Blue';
%             names{i} = int2str(i);
%         end
%         for i = (nseq1+nseq2+1):nseq
%             colors{i} = 'Grey25';
%             names{i} = int2str(i);
%         end
%         adj2pajek(A,names,colors,[],pajekFile);
