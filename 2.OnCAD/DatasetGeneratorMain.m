clear all
save_to_file = true;
% dims = [30,40,50,60,70];

dims = [2, 3, 5, 10, 20];
% dims = [2];
% dims = [3];
% dims = [5];
% dims = [10];
% dims = [20];
% dims = [30];

% dims = [20];
% weights = [0.1, 0.15, 0.2, 0.25, 0.3, 0.5, 1];
weights = [0.1 0.2 0.5 1];

anomaly_rates = [0 0.05];
% anomaly_rates = [0.01];

% insert_noise_after = 100;
max_pop_approximate = 500000;
% num_of_states = 50;


% mean_max = [50000,50000,50000, 50000, 50000];
mean_max = [1500,550,400,380, 400];
% mean_max = [1500];
% mean_max = [1000];
% mean_max = [550];
% mean_max = [400];
% mean_max = [380];
% mean_max = [400];

min_pop_size = 450 + randi(50);
% min_pop_size = 100;

rounds = 1;
with_transition = false;

% iGlobal = 1;
for iGlobal =1:length(dims)
    dim = dims(iGlobal);
    jGlobal = 1;
    num_of_underlying_dists = 100 + randi(10)
    
    
    dist_means = [];
    i = 1;
    k = 1;
    k_counter = 0;
    
    mean_max_temp = mean_max(iGlobal);
    
    covs = cell(1, num_of_underlying_dists);
    while i <= num_of_underlying_dists
        dist_means(i, :) = getRandomMeans(mean_max_temp, floor(mean_max_temp), dim);
        
        covs{i} = generateSPDmatrix(dim)*70;
        %             covs{i} = generateSparseSPDMatrix(dim, 1)*70;
        Stuck = false;
        for j = 1:i-1
            d = FindDistance(dist_means(j,:), inv(covs{i}), dist_means(i,:));
            
            if d < 3*chi2inv(0.9999, dim)
                k = k + 1;
                i = i - 1;
                %                     dim
                %                     k
            end
            if k > 500
                Stuck = true;
                k_counter = k_counter + 1;
                if k_counter > 15
                    if i < num_of_underlying_dists/2
                        mean_max_temp = floor(mean_max_temp + 0.05 * mean_max_temp);    
                    else
                        mean_max_temp = floor(mean_max_temp + 0.02 * mean_max_temp);
                    end
                    'mean_max_temp increased to'
                    mean_max_temp 
                    k_counter = 0;
                end
                
                'distributions not found ...'
                i
                break;
            end
        end
        if Stuck == true
            dist_means = [];
            k = 1;
            i = 0;
        end
        i = i + 1;
    end
    %          ms = cell2mat(dist_means);
    %          plot(dist_means(:,1), dist_means(:,2), '*');
    
    
    state_seq = {};
    %         state_seq = cell(num_of_states, 1);
    mixture_weight = {};
    %         mixture_weight = cell(num_of_states, 1);
    pop_sizes = [];
    %         pop_sizes = zeros(1,num_of_states );
    
    i = 1;
    for jGlobal = 1:length(weights)
        %     while jGlobal <= length(weights)
        min_weight = weights(jGlobal);
        
        while sum(pop_sizes) <= max_pop_approximate
            t_num_of_states = randi(floor(1/min_weight));
            t_state = randsample(num_of_underlying_dists, t_num_of_states)';
            state_seq{i} = t_state;
            weight_remain = 1-t_num_of_states*min_weight;
            mix = ones(1,t_num_of_states) * min_weight;
            mixture_weight{i} = AssignWeight(weight_remain, t_num_of_states, mix);
            pop_sizes(i) = (min_pop_size  + randi(400)*t_num_of_states);
            i = i + 1;
        end
        
        
        %         oneRoundPopSize = sum(pop_sizes);
        %         sumPopSiz = sum(pop_sizes);
        
        
        
        SS = {};
        'generate dataset...'
        for kGlobal = 1:length(anomaly_rates)
            anomaly_rate = anomaly_rates(kGlobal);
            parfor counter = 1:4
                DataGeneratorGeneral(dist_means, covs, pop_sizes, state_seq,mixture_weight , min_weight, dim, anomaly_rate, rounds, with_transition, save_to_file, counter, num_of_underlying_dists);
            end
        end
        jGlobal
    end
    %     iGlobal = iGlobal + 1
end