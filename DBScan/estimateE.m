function [ e ] = estimateE( dataset,  minPts)
%ESTIMATEE Summary of this function goes here
%   Detailed explanation goes here

% dim = size(dataset,2);
n = size(dataset, 1);

dataset= unique(dataset,'rows','stable');
plot_figures = false;

minx = min(dataset(:,1));
maxx = max(dataset(:,1));
miny = min(dataset(:,2));
maxy = max(dataset(:,2));

dists = [];

% Initialization
Points = LOF(dataset, minPts);
adj_mat = zeros(n);
i = 1;
while i <= n
    k = 1;
    for j = Points(i).knn
        adj_mat(i,j) = Points(i).kdist(k);
        dists = [dists , adj_mat(i,j)];
        k = k + 1;
    end
    i = i + 1;
end
if plot_figures
    figure 
    h = histogram(dists);
    vals = h.Values;
    edges = h.BinEdges
else 
   [vals, edges] = histcounts(dists);    
end


diffs = [];


for i = 1:size(vals,2) - 1
    diffs(i) = vals(i) - vals(i+1);
end

e = [];

[val, ind] = max(diffs);
% diffs(ind) = [];
e = edges(ind+1);


end

