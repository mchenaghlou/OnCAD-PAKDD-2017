%======================================================
% Author: Milad Chenaghlou
% Created: 2017-02-02
% Incremental clustering method
%
% Parameters ------------------------------
% dataset = n*dim dataset
% dim = data dimension
% plotFigs = if 0 plots nothin, if 1 plots figures for 2D case.
% filename = file name to store the video
% cluster_similarity_threshold = the threshold for merging clusters, it is 
%   a number in range (0, inf);
% cluster_boundary = the value of P in the Mahalanobis distance formula;
% E = the radius for performing DBScan;
% min_wgh =  the minimum weight of a component in a mixture-distribution 
%   that can be detected by OnCAD;
% remove_cycle = an optional parameter to remove unwanted clusters at these
%   periods
%=======================================================

function [ AnomalyScores clusters ClusterIndexes] = OnCADMain( dataset, dim , plotFigs, filename, cluster_similarity_threshold,cluster_boundary, E, min_wgh, remove_cycle)
% Size of dataset
dataset_size = size(dataset,1);

% Confidence is 1 - alpha
alpha = 0.05;

% Number of observations for each cluster.
means_distance = chi2inv(0.5,dim);
sample_size = MultiVariateNormalSampleSize(alpha, dim, means_distance);

% If a candidate cluster has members above this number, it is confident
% cluster.
tempD = chi2inv(0.1,dim);
confident_cluster_threshold = MultiVariateNormalSampleSize(alpha, dim, tempD);

% Calculate cluster indices
ClusterIndexes = cell(1, dataset_size);
for i = 1:dataset_size
    ClusterIndexes{i} = 0;
end

% Compute buffer_size and initialize buffer
% buffer_size = int16(sample_size * (1/min_wgh)) + 0.1*int16(sample_size * (1/min_wgh));
buffer_size = ceil(sample_size * (1/min_wgh)) + 1;
buffer_size = ceil(buffer_size + buffer_size * 0.05);
active_clusters = [];

%-------------------------------------------

% Initialize cluster vectors
candidate_clusters = {};
ellipse_handles = {};

% Initialize i,j indices
i = sample_size-2; % Index of update cell
j = 0; % Index of check cell

init_cluster = [];

if plotFigs == true
    % anomalyFig = figure;
    writerObj = VideoWriter(filename); % Name it.
    writerObj.FrameRate = 30; % How many frames per second.
    open(writerObj);
end

% The memory associated with the buffer.
A = zeros(1,dataset_size);

% The main loop.
% In this loop, a buffer is simulated with i,j. i is the front and j is
% the back of the buffer.
while i < dataset_size
    if rem(i,1000) == 0
        i
    end
    if i - j + 1 >= buffer_size
        j = j + 1;
    end
    i = i + 1;
    
    % The initialization to form the first cluster.
    if i == sample_size -1
        A(1:sample_size-1) = 1;
        tdataset = dataset(1:sample_size-1,1:dim);
        cov_mat = cov(tdataset );
        inv_cov = inv(cov_mat);
        
        m = mean(tdataset );
        init_cluster.mean =  m;
        init_cluster.chr_mat = inv_cov;
        init_cluster.num_of_members = sample_size - 1;
        init_cluster.last_update = sample_size;
        init_cluster.alpha = sample_size;
        init_cluster.beta = sample_size;
        candidate_clusters = [candidate_clusters, init_cluster];
        
        if E < 0
            x = dsearchn(tdataset,m);
            p2d = pdist2(tdataset, tdataset(x, :));
            p = sort(p2d);
%             E = mean(p(floor(length(p))/2:end))
            E = ceil(max(p2d));
            tempE = E;
            clear temp_dataset seed p2d
        else
            tempE = E;
        end
                


        ClusterIndexes(1:sample_size-1) = {1};
        if plotFigs == true
            hold on
            axis equal
            plot(dataset(1:i,1),dataset(1:i,2),'.');
            
            e = Ellipse_Plot(init_cluster.chr_mat, init_cluster.mean,1,[0,0,0], '-', 2.5);
            ellipse_handles = [ellipse_handles,{e}];
            
        end
        continue;
    end
    
    curr_ob = dataset(i,1:dim);
    if plotFigs == true
        hold on
        if dim == 2
            plot(curr_ob(:,1), curr_ob(:,2), '.');
        end
    end
    % UpdateCell Operations -----------------------------------------------
    
    % 1. Update clusters with weights
    mahal_dists = [];
    member_clusters = [];
    k = 1;
    while k <= size(candidate_clusters,2)
        clus = candidate_clusters{k};
        mahal_dist = FindDistance(curr_ob, clus.chr_mat, clus.mean);
        if mahal_dist < 0
            milad = 1;
        end
        if mahal_dist <= chi2inv(cluster_boundary, dim)
            A(i) = 1;
            member_clusters = [member_clusters, k];
            mahal_dists = [mahal_dists, mahal_dist];
        end
        k = k + 1;
    end
    
    if A(i) == 0
        if plotFigs == true
            hold on
            if dim == 2
                plot(curr_ob(:,1), curr_ob(:,2), '+');
            end
        end
    else
        % Calculate Membership
        alld= 1 ./ (mahal_dists .^ 2);
        weights = alld ./ sum(alld);

        [val, ind] = max(weights);
        ClusterIndexes(i) = {member_clusters(ind)};
        
    end
    
    if length(member_clusters) > 1
        milad = 1;
    end
    
    k = 1;
    while k <= length(member_clusters)
        candidate_clusters{member_clusters(k)} = updateWeightedMean(candidate_clusters{member_clusters(k)}, curr_ob, weights(k));
        candidate_clusters{member_clusters(k)} = updateWeightedChrMat(candidate_clusters{member_clusters(k)}, curr_ob, weights(k));
        
        candidate_clusters{member_clusters(k)}.alpha = int64(candidate_clusters{member_clusters(k)}.alpha) + int64(weights(k));
        candidate_clusters{member_clusters(k)}.beta =  int64(candidate_clusters{member_clusters(k)}.beta) + int64(weights(k).^2);
        
        candidate_clusters{member_clusters(k)}.last_update = i;
        candidate_clusters{member_clusters(k)}.num_of_members = candidate_clusters{member_clusters(k)}.num_of_members + 1;
        
        if plotFigs == true
            if candidate_clusters{member_clusters(k)}.num_of_members > confident_cluster_threshold
                if ishandle(ellipse_handles{member_clusters(k)}{2}) delete(ellipse_handles{member_clusters(k)}{2}); end
                if ishandle(ellipse_handles{member_clusters(k)}{1}) delete(ellipse_handles{member_clusters(k)}{1}); end
                e = Ellipse_Plot(candidate_clusters{member_clusters(k)}.chr_mat, candidate_clusters{member_clusters(k)}.mean,1,[0,0,0], '-', 2.5);
                ellipse_handles(member_clusters(k)) = {e};
            else
                if (ishandle(ellipse_handles{member_clusters(k)}{2})) delete(ellipse_handles{member_clusters(k)}{2}); end
                if (ishandle(ellipse_handles{member_clusters(k)}{1})) delete(ellipse_handles{member_clusters(k)}{1}); end
                e = Ellipse_Plot(candidate_clusters{member_clusters(k)}.chr_mat, candidate_clusters{member_clusters(k)}.mean,1,[0,0,1], '-.', 2);
                ellipse_handles(member_clusters(k)) = {e};
            end
        end
        
        k = k + 1;
    end
    
    
%     % 2. Remove clusters that do not get updated recently.
%     cc_counter = 1;
%     %     active_clusters = candidate_clusters;
%     if rem(i, 50) == 0
%         while cc_counter <= length(active_clusters)
%             try
%                 candidate_cluster = candidate_clusters{active_clusters(cc_counter)};
%             catch ME
%                 msgText = getReport(ME)
%                 error('milad')
%             end
%             if candidate_cluster.last_update < j && candidate_cluster.num_of_members < confident_cluster_threshold
%                 try
%                     ClusterIndexes(candidate_clusters{active_clusters(cc_counter)}.indexes) = {0};
%                 catch ME
%                     msgText = getReport(ME)
%                     error('milad')
%                 end
%                 
%                 candidate_clusters(active_clusters(cc_counter)) = [];
%                 pivot = active_clusters(cc_counter);
%                 active_clusters(cc_counter) = [];
%                 for k = cc_counter:length(active_clusters)
%                     if active_clusters(k) > pivot
%                         active_clusters(k) = active_clusters(k) - 1;
%                     end
%                 end
%                 
%                 if plotFigs == true
%                     if (ishandle(ellipse_handles{cc_counter}{2})) delete(ellipse_handles{cc_counter}{2}); end;
%                     if (ishandle(ellipse_handles{cc_counter}{1})) delete(ellipse_handles{cc_counter}{1}); end;
%                     ellipse_handles(cc_counter) = [];
%                 end
%                 cc_counter = cc_counter -1;
%                 A(j) = 0;
%                 
%             end
%             cc_counter = cc_counter + 1;
%         end
%     end
    
    
    
    
%     3. Merge Clusters - 2 cases : Candidate-clusters within each other are
%     merged into 1. Candidate and confident clusters that are very similar
%     get removed.
%     if rem(i, remove_cycle) == 0
%         k = 1;
%         while k <= size(candidate_clusters,2) - 1
%             kk = k + 1;
%             cur_clu = candidate_clusters{k};
%             while kk <= size(candidate_clusters,2)
%                 next_clus = candidate_clusters{kk};
%                 d = FocalDistancemd(cur_clu.chr_mat, cur_clu.mean, next_clus .chr_mat, next_clus .mean );
%                 if d < cluster_similarity_threshold
%                     if plotFigs == true
%                         if (ishandle(ellipse_handles{kk}{2})) delete(ellipse_handles{kk}{2}); end;
%                         if (ishandle(ellipse_handles{kk}{1})) delete(ellipse_handles{kk}{1}); end;
%                         ellipse_handles(kk) = [];
%                     end
%                     candidate_clusters(kk) = [];
%                     kk = kk - 1;
%                 end
%                 kk = kk + 1;
%             end
%             k = k + 1;
%         end
%     end

    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CheckCell Operations
    % If the observation is not assigned to any state-tracker, do
    % clustering over the buffer.
    % If new clusters are formed, create state-trackers for each
    % else : report anomaly.
    
    if j > 0
        if A(j) == 0
            % 1. Do clustering over observations not assigned to any state-tracker
            free_inds = find(A(j:i) == 0) + j - 1;
            if length(free_inds) >= sample_size -2
                frees = dataset(free_inds, 1:dim);
                %%% ------------Clustering algorithm --------------------------------------
                %%%%%%%%%%%%%%%%%%%%%%%%%%% DB SCAN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                [dbClu, ptsC, centres] = dbscan(frees', tempE , sample_size );
%                 if(size(dbClu, 2) == 0)
%                     tempE = tempE + tempE*0.05;
%                 else
%                     tempE = E;
%                 end
                clu = {};
                
%             if size(dbClu,2) == 0
% %                  i
% %                  s = sprintf('DBScan could not find clusters with parameters: E = %d, MinPts = %d', E, sample_size);
% %                  s
% %                  seed = frees(1,:);
% %                  p2d = pdist2(frees, seed);
% %                  s = sortrows(p2d);
%             end
                for counter = 1:size(dbClu,2)
                    c.mean = centres(:,counter)';
                    c.indexes = free_inds(dbClu{counter});
                    if length(c.indexes) >= sample_size
                        cc = cov(dataset(c.indexes, 1:dim));
%                         cc = [inf, inf; inf, inf]
                        if sum(sum(isinf(cc))) > 0 || sum(sum(isnan(cc))) > 0 || sum(sum(abs(cc) < 10^(-15))) ...
                                || sum(sum((cc > 0))) == 0
                            milad = 1;
                        else
                            c.chr_mat = inv(cc);
%                             c.chr_mat = cc;
                            c.last_update = max(c.indexes);
                            A(c.indexes) = 1;
                            ClusterIndexes(c.indexes) = {size(candidate_clusters,2)+1};

                            clu = [clu, c];                            
                        end

                    end
                    
                end
                
                %%% -----------------------------------------------------------
                % use weighted batch DCAD model to create chr mat.
                for k = 1:size(clu,2)
                    clust = clu{k};
                    clust.num_of_members = length(clust.indexes);
                    clust.beta = length(clust.indexes);
                    clust.alpha = length(clust.indexes);
                    candidate_clusters = [candidate_clusters, clust];
%                     active_clusters  = [active_clusters , length(candidate_clusters)];
                    if plotFigs == true
                        hold on
                        e = Ellipse_Plot(clust.chr_mat, clust.mean,1,[0,1,0], '-', 2.5);
                        ellipse_handles = [ellipse_handles, {e}];
                    end
                    
                end
            else
                % What do we do here?
                % Nothing.
                
            end
        end
        
    end

    
    if plotFigs == true
        %% Set up the movie.
        %     figure(normal); % Makes sure you use your desired frame.
        if mod(i,10)==0 % Uncomment to take 1 out of every 4 frames.
            frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
            writeVideo(writerObj, frame);
        end
        
        %end
    end
    
    if isempty(ClusterIndexes{i})
        error('something is wrong.')
    end
    %          elapse = toc;
    %          times = [times, elapse];
end
if plotFigs == true
    close(writerObj); % Saves the movie.
end
AnomalyScores = 0;
clusters = candidate_clusters;
end


