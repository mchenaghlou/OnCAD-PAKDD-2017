function [ curr_cluster ] = updateWeightedChrMat( curr_cluster, curr_ob, w)
%UPDATEWEIGHTEDCHRMAT Summary of this function goes here
%   Detailed explanation goes here
% weights parameter is to be implemented


% n = curr_pot_cluster.num_of_members + 1;
% dim = size(curr_ob,2);
% sn_1 = curr_pot_cluster.chr_mat;
% mn_1 = curr_pot_cluster.mean';
% x_n = curr_ob';
% 
% 
% enumerator = (((x_n - mn_1)*(x_n - mn_1)' * sn_1 ));
% denominator = ((n*n - 1)/n) + ((x_n - mn_1)' * sn_1) * (x_n - mn_1);
% coef = ((n * sn_1) / (n-1));
% chr_mat = coef * (diag(ones(1,dim))  - enumerator/ denominator);

curr_ob = curr_ob';
mK = curr_cluster.mean';

sK_1 = curr_cluster.chr_mat;

betaK = int64(curr_cluster.beta);
betaKplus1 = int64(curr_cluster.beta) + w.^2;
alphaK = int64(curr_cluster.alpha);
alphaKplus1 = (curr_cluster.alpha) + w;
if betaK > 88
    milad = 1;
    a = 88*7832;
end
b = (int64(betaKplus1.^2) - int64(alphaKplus1));
khi_enum = int64(betaK) .* int64(b);
khi_denom = (int64(betaKplus1)*(int64(betaK.^2) - int64(alphaK)));
khiK = double(khi_enum) / double(khi_denom);

delta_enum  = int64(betaKplus1) * int64((betaK .^2 - alphaK));
delta_denom = int64(betaK) * w * int64((betaKplus1 + w -2));
delta = double(delta_enum) / double(delta_denom);

enumerator = ((sK_1 * (curr_ob - mK)) * ((curr_ob - mK)' * sK_1));
denominator = (delta + ((curr_ob - mK)' * sK_1 * (curr_ob - mK)));
new_sn_1 = khiK * (sK_1  - enumerator/denominator);
curr_cluster.chr_mat = new_sn_1;
end

