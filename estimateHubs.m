function k = estimateHubs(TransNet)
% AM = symmetrize(TransNet) > 0;
AM = (TransNet + TransNet') > 0;
n = size(AM,1);
% deg = degrees(AM);
deg = sum(AM);
deg = sort(deg,'descend');

T = clusterdata(deg','linkage','average','maxclust',2);
c = T(1);
k = size(find(T == c),1);

% pmin = 100;
% mmax = 0;
% k = 0;
% for  i = 2:n
% %     group = zeros(1,n);
% %     group(i:end) = 1;
% %     p = anova1(deg,group,'off');
% %     if p < pmin
% %         pmin = p;
% %         k = i-1;
% %     end
%     m1 = mean(deg(1:(i-1)));
%     m2 = mean(deg(i:end));
%     m = m1/m2;
%     if m > mmax
%         mmax = m;
%         k = i-1;
%     end
% end

% Esmax = buildSmaxGraph(degrees(AM));
% AMsmax= edgeL2adj(Esmax);
% AMsmax = symmetrize(AMsmax);
% sm = sMetric(AMsmax)/2;
% sm1 = sMetric(AM)/2;
% names_pat = cellfun(@num2str, num2cell(1:n), 'UniformOutput', false);
% bg = biograph(sparse(AMsmax),names_pat,'ShowWeights','on');
% view(bg);

k

