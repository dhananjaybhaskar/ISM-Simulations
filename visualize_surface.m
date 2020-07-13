function [] = visualize_surface(X, itr, vis_x, vis_y, vis_z, x_limits, y_limits, z_limits)

    fig = figure('visible', 'off');
    mesh(vis_x, vis_y, vis_z, 'edgecolor', 'k');
    alpha 0.75;
    daspect([1 1 1])
    xlim(x_limits)
    ylim(y_limits)
    zlim(z_limits)
    hold on;
    scatter3(X(:,1), X(:,2), X(:,3));
    plot3(path_matrix_x, path_matrix_y, path_matrix_z, 'r', 'LineWidth', 2);
    fname = strcat('sim_', sprintf('%03d',itr), '.png');
    saveas(fig, fname, 'png');
    
end