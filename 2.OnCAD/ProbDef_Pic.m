clear all
hold on
axis equal

a1 = 0;
b1 = 14;
a2 = 1;
b2 = 10;

axis([a1 b1 a2 b2])
mu1 = [2,3];
sigma1 = [0.3,0.05;0.05,0.3];
r = mvnrnd(mu1,sigma1,1000);
dataset = r;

mu2 = [5,6];
sigma2 = [0.3,0.1;0.1,0.3];
r = [r; mvnrnd(mu2,sigma2,500);]
dataset = [dataset;r];

mu3 = [5,3];
sigma3 = [0.3,0.07;0.07,0.3];
r = [r; mvnrnd(mu3,sigma3,200);]
dataset = [dataset;r];
plot(r(:, 1), r(:, 2), '.')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu4 = [10,3];
sigma4 = [0.3,0.1;0.1,0.3];
r = mvnrnd(mu4,sigma4,300);
dataset = [dataset;r];

mu5 = [12,5.5];
sigma5 = [0.3,0.1;0.1,0.3];
r = [r; mvnrnd(mu5,sigma5,1000);]
dataset = [dataset;r];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

anoms = rand(100,2);

anoms(:, 1) = (b1-a1) * anoms(:, 1) + a1;

anoms(:, 2) = (b2-a2) * anoms(:, 2) + a2;

r = [r; anoms]
plot(r(:, 1), r(:, 2), '.')
[AnomalyScores, Clusters, ClusterIndexes] = OnCAD( ...
    dataset, 2, false, 1, 0.99, -1, 0.1, 100000);

for m = 1:length(Clusters)
    clust = Clusters{m};
    Ellipse_Plot(clust.chr_mat, clust.mean,1,[0,0,0], '-', 2.5);
end




mu6 = [9,8];
sigma6 = [0.3,0.1;0.1,0.3];
r = mvnrnd(mu6,sigma6,100);
dataset = [dataset;r];

plot(r(:, 1), r(:, 2), '.')
