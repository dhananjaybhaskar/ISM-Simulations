function visualize_trajectories(X, itr, previous_steps, colors, vis_x, vis_y, vis_z, x_limits, y_limits, z_limits)

    fig = figure('visible', 'off'); 
    
    h = mesh(vis_x, vis_y, vis_z, 'EdgeColor', 'none', 'FaceColor', [0.7, 0.7, 0.7], 'FaceAlpha', 1, 'linestyle', '-');
    camlight
    h.FaceLighting = 'gouraud';
    h.AmbientStrength = 0.3;
    h.DiffuseStrength = 0.8;
    h.SpecularStrength = 0.9;
    h.SpecularExponent = 25;
    h.BackFaceLighting = 'lit';
    
    alpha 1
    daspect([1 1 1])
    xlim(x_limits)
    ylim(y_limits)
    zlim(z_limits)
    hold on;
    
    scatter3(X(:,1), X(:,2), X(:,3), 30, 'filled');
    for i = 1:length(X)
        plot3(previous_steps(i,:,1), previous_steps(i,:,2), previous_steps(i,:,3), 'LineWidth', 1.2, 'Color', colors(i, :));
    end
    
    zoom(1.2)
    set(gca, 'visible', 'off')
    fname = strcat('sim_', sprintf('%03d', itr), '.png');
    saveas(fig, fname, 'png');
    close
    
end