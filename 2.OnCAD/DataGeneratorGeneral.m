function [ S ] = DataGeneratorGeneral(states, covs, pop_sizes, states_seq, mixture_weights, min_weight, dim, anomaly_rate, rounds, with_transition, save_to_file, iter, num_of_underlying_dists)
%DATAGENERATOR3 Summary of this function goes here
%   Detailed explanation goes here
% range default



    
    fileName = sprintf('C:/Users/mchenaghlou/Documents/MyPapers/2.IncClustering/Data/DataSets-Synthetic-4/dataset-dim%s-weight%s-anomalyRate%s-iter%d.csv',num2str(dim), num2str(min_weight), num2str(0), iter);
    if exist(fileName, 'file') == 2
        S = csvread(fileName);    
        fileName = sprintf('C:/Users/mchenaghlou/Documents/MyPapers/2.IncClustering/Data/DataSets-Synthetic-4/label-dim%s-weight%s-anomalyRate%s-iter%d.csv',num2str(dim), num2str(min_weight), num2str(0), iter);
        label = csvread(fileName);
        S = insertAnomalies(S, label, anomaly_rate, dim);
    else
        S = [];
        for round_counter = 1:rounds
            S_temp = generateCleanDataset(states, covs, pop_sizes, states_seq, mixture_weights, min_weight, dim, anomaly_rate, rounds, with_transition, save_to_file, iter, num_of_underlying_dists);    
            S = [S; S_temp];
        end
        
        S = insertAnomalies(S(:, 1:dim), S(:, dim+1), anomaly_rate, dim);
    end
    
    

    
 
    % minX = min(S(:,1));
    % maxX = max(S(:,2));
    % minY = min(S(:,1));
    % maxY = max(S(:,2));
    % buff = 50;
    %
    %
    % writerObj = VideoWriter('dataset-2d.avi'); % Name it.
    % writerObj.FrameRate = 20; % How many frames per second.
    % open(writerObj);
    % hold on
    % for i = 1:floor(wholePopSize/buff )
    % %     plot3(S((i-1)*10 + 1:i*10,1),S((i-1)*10 + 1:i*10, 2), S((i-1)*10 + 1:i*10 ,3), '.');
    %     chunk = S((i-1)*buff  + 1:i*buff ,:);
    %     anoms = chunk(chunk(:, 3) == 0, :);
    %     norms = chunk(chunk(:, 3) ~= 0, :);
    %     hold on
    %     plot(norms(:,1),norms(:,2), '.g');
    %     plot(anoms(:,1),anoms(:,2), '+r');
    % %     axis([minX maxX minY maxY])
    %     pause(0.01);
    %
    % %     frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
    % %     writeVideo(writerObj, frame);
    %
    % end
    % close(writerObj); % Saves the movie.
    
    % save('S', 'S');
    



plot(S(:,1), S(:,2), '.');
axis equal
if(size(S,1) < 500000)
    'error ..............'
end
if save_to_file == true
    toSave = sprintf('./../Data/DataSets-Synthetic-4/summary-dim=%d-weight=%s-anomalyRate=%s-iter=%d.mat', dim, num2str(min_weight), num2str(anomaly_rate), iter);
    fileName=sprintf('./../Data/DataSets-Synthetic-4/dataset-dim%s-weight%s-anomalyRate%s-iter%s.csv', num2str(dim), num2str(min_weight), num2str(anomaly_rate), num2str(iter));
    labelName=sprintf('./../Data/DataSets-Synthetic-4/label-dim%s-weight%s-anomalyRate%s-iter%s.csv', num2str(dim), num2str(min_weight), num2str(anomaly_rate), num2str(iter));
    save(toSave ,'states', 'covs', 'pop_sizes', 'states_seq', 'mixture_weights', 'dim', 'anomaly_rate', 'rounds', 'with_transition', 'num_of_underlying_dists');
    %     save(fileName ,'S')
    
    
    csvwrite(fileName, S(:, 1:dim));
    csvwrite(labelName, S(:, dim+1));
end




clear alpha cc co iii inds init iter k ls maxK min_weight sample_size save_to_file t tcurcov tcurmean toSave tS anomalies after before m mas mis sp toInsert anomaly_rate ans covs cur_cov cur_dist cur_mean cur_mixture cur_pop_size curr_dists dim i i_temp ii indices j_temp ii indices j_temp mixture_weights move_states ms n nomin numOfAnomalies obs  oneRoundPopSize path path_length pop_sizes pre_dists rounds states states_seq steps sum_path_lengths this_transit_length transit_dists transit_dest transit_length transit_source wholePopSize with_transition


end

