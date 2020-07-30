function [particle_indices] = find_neighbors(X, mesh_idxes, mesh_dists, dist_mat, source_pt_idx, threshold)

    particle_indices = zeros(length(X), 1);
    source_dist = mesh_dists(source_pt_idx);
    source_mesh_idx = mesh_idxes(source_pt_idx);
    for i = 1:length(X)
        if(i ~= source_pt_idx)
            point_dist = mesh_dists(i);
            point_mesh_idx = mesh_idxes(i);
            total_dist = source_dist + dist_mat(source_mesh_idx, point_mesh_idx) + point_dist;
            if(total_dist <= threshold)
                particle_indices(i) = i;
            end
        end
    end
    
    particle_indices = particle_indices(particle_indices ~= 0);
    
end