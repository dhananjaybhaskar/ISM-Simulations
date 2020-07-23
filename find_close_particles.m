function [particle_indices] = find_close_particles(X, mesh_x, mesh_y, mesh_z, dist_mat, source_pt_idx, threshold)
    [base_dist, source_mesh_idx] = closest_mesh_point(X(source_pt_idx, :), mesh_x, mesh_y, mesh_z);
    particle_indices = zeros(numel(X),1);
    for i = 1:numel(X)
        if(i ~= source_pt_idx)
            [point_dist, point_mesh_idx] = closest_mesh_point(X(i, :), mesh_x, mesh_y, mesh_z);
            total_dist = base_dist + dist_mat(source_mesh_idx, point_mesh_idx) + point_dist;
            if(total_dist <= threshold)
                particle_indices(i) = i;
            end
        end
    end
    particle_indices = particle_indices(particle_indices ~= 0);
end