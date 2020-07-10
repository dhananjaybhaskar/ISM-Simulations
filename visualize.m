function [] = visualize(X, itr, pt_1_idx, pt_2_idx,x,y,z,num_cols, next, x_limits, y_limits, z_limits)
    point_1 = X(pt_1_idx, :);
    point_2 = X(pt_2_idx, :);
    [point_1_dist, point_1_mesh_idx] = closest_mesh_point(point_1,x,y,z);

    [point_2_dist, point_2_mesh_idx] = closest_mesh_point(point_2,x,y,z);
    

    translated_path_idxes = compute_path(next,point_1_mesh_idx,point_2_mesh_idx,num_cols);
    len = length(translated_path_idxes);
    path_matrix_x = zeros(len + 2, 1);
    path_matrix_y = path_matrix_x;
    path_matrix_z = path_matrix_x;
    
    path_matrix_x(1) = point_1(1);
    path_matrix_y(1) = point_1(2);
    path_matrix_z(1) = point_1(3);
    
    for i = 2:(length(translated_path_idxes) + 1)
        path_matrix_x(i) = x(translated_path_idxes(i - 1, 1), translated_path_idxes(i - 1, 2));
        path_matrix_y(i) = y(translated_path_idxes(i - 1, 1), translated_path_idxes(i - 1, 2));
        path_matrix_z(i) = z(translated_path_idxes(i - 1, 1), translated_path_idxes(i - 1, 2));
    end
    
    path_matrix_x(len + 2) = point_2(1);
    path_matrix_y(len + 2) = point_2(2);
    path_matrix_z(len + 2) = point_2(3);

    fig = figure('visible', 'off');
    mesh(x, y, z, 'edgecolor', 'k');
    alpha 0.75;
    daspect([1 1 1])
    xlim(x_limits)
    ylim(y_limits)
    zlim(z_limits)
    hold on;
    scatter3(X(:,1), X(:,2), X(:,3));
    plot3(path_matrix_x, path_matrix_y, path_matrix_z, 'r', 'LineWidth',2);
    fname = strcat('sim_', sprintf('%03d',itr), '.png');
    saveas(fig, fname, 'png');

end