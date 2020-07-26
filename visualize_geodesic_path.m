function [] = visualize_geodesic_path(X, itr, pt_1_idxs, pt_2_idxs, vis_x, vis_y, vis_z, mesh_x, mesh_y, mesh_z, ...
                                        num_cols, next, x_limits, y_limits, z_limits)

    num_points = numel(pt_1_idxs);
    assert(num_points == numel(pt_2_idxs), "Array size does not match.");
    assert(num_points >= 1, "Please specify atleast one pair of points.");
        
    fig = figure('visible', 'off');
    mesh(vis_x, vis_y, vis_z, 'edgecolor', 'k');
    alpha 0.75;
    daspect([1 1 1])
    xlim(x_limits)
    ylim(y_limits)
    zlim(z_limits)
    hold on;
    scatter3(X(:,1), X(:,2), X(:,3), 30, 'bo');
    
    for pti = 1 : num_points
                                    
        point_1 = X(pt_1_idxs(pti), :);
        point_2 = X(pt_2_idxs(pti), :);

        [~, point_1_mesh_idx] = closest_mesh_point(point_1, mesh_x, mesh_y, mesh_z);
        [~, point_2_mesh_idx] = closest_mesh_point(point_2, mesh_x, mesh_y, mesh_z);

        translated_path_idxes = compute_path(next, point_1_mesh_idx, point_2_mesh_idx, num_cols);
        len = length(translated_path_idxes);

        path_matrix_x = zeros(len + 2, 1);
        path_matrix_y = path_matrix_x;
        path_matrix_z = path_matrix_x;

        path_matrix_x(1) = point_1(1);
        path_matrix_y(1) = point_1(2);
        path_matrix_z(1) = point_1(3);

        for i = 2:(length(translated_path_idxes) + 1)
            path_matrix_x(i) = mesh_x(translated_path_idxes(i - 1, 1), translated_path_idxes(i - 1, 2));
            path_matrix_y(i) = mesh_y(translated_path_idxes(i - 1, 1), translated_path_idxes(i - 1, 2));
            path_matrix_z(i) = mesh_z(translated_path_idxes(i - 1, 1), translated_path_idxes(i - 1, 2));
        end

        path_matrix_x(len + 2) = point_2(1);
        path_matrix_y(len + 2) = point_2(2);
        path_matrix_z(len + 2) = point_2(3);

        plot3(path_matrix_x, path_matrix_y, path_matrix_z, 'r', 'LineWidth', 1.0);
        scatter3(point_1(1), point_1(2), point_1(3), 30, 'ro', 'filled')
        scatter3(point_2(1), point_2(2), point_2(3), 30, 'ro', 'filled')
        
    end
    
    hold off;
    fname = strcat('sim_', sprintf('%03d', itr), '.png');
    saveas(fig, fname, 'png');
    close
    
end