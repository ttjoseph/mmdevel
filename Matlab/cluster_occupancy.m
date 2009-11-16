function [occupancy] = cluster_occupancy(clusters, halfwindowsize)
% Use a sliding window for graphing cluster occupancy
% That is, how many cluster IDs are mentioned in the window
num_atoms = size(clusters, 1);
occupancy = zeros(num_atoms, 1);

for i = halfwindowsize + 1 : num_atoms - halfwindowsize
    clusters_seen = zeros(num_atoms, 1);
    
    for w = i - halfwindowsize : i + halfwindowsize - 1
        clusters_seen(clusters(w)) = 1;
    end
    occupancy(i) = sum(clusters_seen);
    %occupancy(i) = std(clusters_seen);

end