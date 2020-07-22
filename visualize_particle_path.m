function visualize_particle_path(X, itr, previous_steps, colors,  vis_x, vis_y, vis_z, x_limits, y_limits, z_limits)

    fig = figure('visible', 'off'); clf;
    mesh(vis_x, vis_y, vis_z, 'edgecolor', 'k');
    alpha 0.75;
    daspect([1 1 1])
    xlim(x_limits)
    ylim(y_limits)
    zlim(z_limits)
    hold on;
    scatter3(X(:,1), X(:,2), X(:,3));
    for i = 1:length(X)
        plot3(previous_steps(i,:,1), previous_steps(i,:,2), previous_steps(i,:,3),'LineWidth', 1.5, 'Color',colors(i, :));
    end
    fname = strcat('sim_', sprintf('%03d',itr), '.png');
    saveas(fig, fname, 'png');
end