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
    fname = strcat('sim_', sprintf('%03d',itr), '.png');
    saveas(fig, fname, 'png');
    
end