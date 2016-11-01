function val = objTransNetQuadDeg2(DM,nhubs,aux)

k = nhubs;
% k = 1;
sigma = 0.5;
epsilon = 0.01;

nvert = size(DM,2);
AM = (DM+DM' > 0);
degseq = sum(AM,2);

% AM = (DM > 0);
% degseq = sum(AM,2);

val = 0;
for i=1:nvert
    for j=(i+1):nvert
        if AM(i,j) == 1
            val = val + degseq(i)*degseq(j);
        end
    end
end

% Esmax = buildSmaxGraph(degrees(AM));
% AMsmax= edgeL2adj(Esmax);
% AMsmax = symmetrize(AMsmax);
% sm = sMetric(AMsmax)/2;
% val = val/sm;

s = sCore(nvert,k);
% s = 171;
val = exp(-abs(val-s)/s);
% val = val/s;
% val = 1 - abs(val-s)/s;
if ~isempty(aux)
    val = (val - aux(1))/(3*aux(2));
    val = (val + 1)/2;
end

% val = val/(nvert*(nvert-1)/2);
% s = s/(nvert*(nvert-1)/2);
% if val <= s
%     val = exp((val-s)/s);
% else
%     val = 0;
% end


% val = 10*abs(val-s)/s;
% p = normcdf([val-epsilon val+epsilon]);
% val = p(2)-p(1);

% possK = 1:(2*k-1);
% l = 3;
% step = 2*l/(2*k-1);
% inter = -l:step:l;
% 
% for i = 1:(2*k-1)
%     
% end



