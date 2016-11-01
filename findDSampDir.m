function DSamp_dir = findDSampDir(DSamp, thr)
nSamp = size(DSamp,1);
AMSamp = zeros(nSamp,nSamp);

for u=1:nSamp
    for v=(u+1):nSamp
            m = min(DSamp(u,v),DSamp(v,u));
            if m<=thr
                if DSamp(u,v) == m
                    AMSamp(u,v) = 1;
                end
                 if DSamp(v,u) == m
                    AMSamp(v,u) = 1;
                end
            end
    end
end

DSamp_dir = DSamp.*AMSamp;