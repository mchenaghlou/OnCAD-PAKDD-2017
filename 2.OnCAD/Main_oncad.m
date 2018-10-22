clear all
load('./2.OnCAD/evolving_S1.mat')
dim = size(dataset, 2);
plotFigs = 1;
filename = 'out8.avi';
min_wgh = 0.1;
[AnomalyScores, Clusters, ClusterIndexes] = OnCAD( dataset, dim, plotFigs, filename, 3,0.99, -1, min_wgh, 500);