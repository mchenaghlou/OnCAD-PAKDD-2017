function [ toRet ] = getRandomMeans( mean_max_1, mean_max_2, dim )
%GETRANDOMMEANS Summary of this function goes here
%   Detailed explanation goes here
    
    toRet = [];
    toRet = [toRet, randi(mean_max_1 ,1 , floor(dim/2))+500];
    toRet = [toRet, randi(mean_max_2 ,1 , dim - floor(dim/2))+500];
    
   
    

end

